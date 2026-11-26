"""ACL's block at ~125 significant digits (python-flint / Arb), so every rung of the ladder can be
resolved down to the bottom (expected ~1e-54 at lam^2 = 12), not just the top few above the
double-precision floor.

Same construction as fourier_block.py:
  Weil's form on ACL's Fourier space E_N of the log window, from 1-D transforms of H
  (notes/autocorrelation.md); conditioning vectors = E-images of ACL's prolate combinations (Fact 5.3),
  projected on E_N; Schur complement on V = span{p_k} + pole.
Everything numeric is done in Arb balls at PREC bits; quadrature is piecewise 48-point Gauss-Legendre
(pieces <= 1.5 periods of the top frequency; geometric grading at U = 0, where H has a
(1 - y) log(1 - y) term), so the quadrature error is far below the quantities we resolve."""
import math, sys, time, json
import numpy as np
from flint import arb, acb, arb_mat, ctx
import mpmath
from mpmath import mp
from mpmath.calculus.quadrature import GaussLegendre

import os
PREC = int(os.environ.get('HP_PREC', '420'))
PIECE = float(os.environ.get('HP_PIECE', '1.5'))     # pieces <= PIECE periods of the top frequency
GRADE = int(os.environ.get('HP_GRADE', '180'))
DIGITS = int(os.environ.get('HP_DIGITS', '0'))      # if set: choose points/grading for ~10^-DIGITS quadrature error
ctx.prec = PREC
mp.prec = PREC
DIG = int(PREC * 0.30103) + 5
ZERO, ONE, HALF = arb(0), arb(1), arb("0.5")
PI = arb.pi()
EG = arb.const_euler()
LOG2PI = (2 * PI).log()


def a_of(x):                      # mpmath -> arb
    return arb(mpmath.nstr(x, DIG, strip_zeros=False))


_GL = {}
def gl(npts):
    """Gauss-Legendre nodes/weights on [-1, 1] as arbs; npts in (12, 24, 48)."""
    if npts not in _GL:
        deg = {12: 3, 24: 4, 48: 5, 96: 6, 192: 7}[npts]
        nodes = GaussLegendre(mp).calc_nodes(deg, mp.prec)
        _GL[npts] = [(a_of(x), a_of(w)) for x, w in nodes]
    return _GL[npts]


def prime_powers(X):
    out = []
    for p in range(2, int(X) + 1):
        if all(p % d for d in range(2, int(p ** 0.5) + 1)):
            q = p
            while q <= X:
                out.append((q, p)); q *= p
    return sorted(out)


class Form:
    def __init__(self, lam2):
        self.lam2 = lam2
        self.L = arb(lam2).log()
        self.pp = prime_powers(lam2)
        self.Lam = {q: arb(p).log() for q, p in self.pp}

    def H(self, U, exclude=()):
        x = U.exp()
        xf = float(x.mid())
        S = ZERO
        for q, p in self.pp:
            if q <= xf and q not in exclude:
                S += self.Lam[q] * (x / q - 1)
        G = S - x * U + (1 + EG) * x
        y = 1 / x
        omy = -((-U).expm1())                         # 1 - y, accurate for small U
        if U.is_zero():
            T = ONE - 2 * arb(2).log() / 2
        else:
            T = ONE - (1 + y) / (2 * y) * y.log1p() + omy / (2 * y) * omy.log()
        return (-U / 2).exp() * (G - LOG2PI + T)

    def ramp(self, U, q):
        x = U.exp()
        v = x / q - 1
        return (-U / 2).exp() * v if float(v.mid()) > 0 else ZERO


def quad_nodes(L, breaks, wmax, grade=True):
    """breaks: arbs in (0, L) (exact kink locations). Pieces <= 1.5 periods of wmax, 48-point GL;
    the first piece is replaced by a geometric grading toward 0 if grade."""
    Lf = float(L.mid())
    pts = [ZERO] + sorted([b for b in breaks if 0 < float(b.mid()) < Lf], key=lambda b: float(b.mid())) + [L]
    hmax = PIECE * 2 * math.pi / max(wmax, 1.0)
    pieces = []
    mainpts = 96 if DIGITS > 100 else 48          # GL error ~ (e w'/4n)^{2n}, w' <= 1.5 pi: 96 points -> 1e-284
    for A, B in zip(pts[:-1], pts[1:]):
        m = max(1, int(math.ceil((float(B.mid()) - float(A.mid())) / hmax)))
        for i in range(m):
            pieces.append((A + (B - A) * i / m, A + (B - A) * (i + 1) / m, mainpts))
    if grade:
        # main pieces close to the singularity at 0 (distance < 20 piece lengths): 96 points (192 at high DIGITS)
        near = 192 if DIGITS > 140 else 96
        pieces = [(A, B, near) if i < 20 else (A, B, n) for i, (A, B, n) in enumerate(pieces)]
        B1 = pieces[0][1]
        g = []
        if DIGITS:
            # ratio-2 pieces, Bernstein rho = 5.83: error ~ 5.83^{-2n} 2^{-i}; n minimal for 10^{-DIGITS-5}
            for i in range(int(DIGITS / 0.301) + 20):
                for npts in (24, 48, 96, 192):
                    if 1.531 * npts + 0.301 * i > DIGITS + 5:
                        break
                g.append((B1 / arb(2) ** (i + 1), B1 / arb(2) ** i, npts))
        else:
            for i in range(GRADE):
                npts = 96 if i < 60 else (48 if i < GRADE // 2 else 24)    # ratio-2 pieces: Bernstein rho = 5.83
                g.append((B1 / arb(2) ** (i + 1), B1 / arb(2) ** i, npts))
        pieces = g + pieces[1:]
    out = []
    for A, B, npts in pieces:
        hw, c = (B - A) / 2, (B + A) / 2
        for xi, wi in gl(npts):
            out.append((hw * xi + c, hw * wi))
    return out


def transforms(nodes, fvals, L, N):
    """S_m = int f sin(w_m s), C2_m = 2 int f (L - s) cos(w_m s), m = 0..N, w_m = 2 pi m / L."""
    S = [ZERO] * (N + 1); C2 = [ZERO] * (N + 1)
    w1 = 2 * PI / L
    for (s, w), f in zip(nodes, fvals):
        fw = f * w; fwl = fw * (L - s) * 2
        s1, c1 = (w1 * s).sin_cos()
        sm, cm = ZERO, ONE
        for m in range(N + 1):
            if m % 32 == 0 and m:                    # restart: the ball recurrence loses up to half a bit per step (wrapping),
                sm, cm = (w1 * s * m).sin_cos()      # which silently destroyed the form above N of about 750 at 420 bits
            S[m] += fw * sm; C2[m] += fwl * cm
            sm, cm = sm * c1 + cm * s1, cm * c1 - sm * s1
    return S, C2


def gram(S, C2, f0, fL, L, N):
    """Even (N+1)^2 and odd N^2 Gram matrices of the form whose kernel has transforms S, C2."""
    sqL = L.sqrt()
    w = [2 * PI * m / L for m in range(N + 1)]
    def Sf(m): return S[m] if m >= 0 else -S[-m]
    def wf(m): return w[m] if m >= 0 else -w[-m]
    const = 2 * (f0 - fL) / L
    cache = {}
    def Qc(m, k):
        key = (m, k)
        if key in cache: return cache[key]
        am = acb(HALF, wf(m)) / sqL; ak = acb(HALF, wf(k)) / sqL
        Jm, Jk = acb(0, 2 * Sf(m)), acb(0, 2 * Sf(k))
        I = C2[abs(m)] if m == k else 2 * (Sf(k) - Sf(m)) / (wf(m) - wf(k))
        v = am * ak.conjugate() * I + (am * Jm + (ak * Jk).conjugate()) / sqL + const
        cache[key] = v
        return v
    Qe = arb_mat(N + 1, N + 1); Qo = arb_mat(N, N)
    r2 = arb(2).sqrt()
    for a in range(N + 1):
        for b in range(a, N + 1):
            sg = 1 if (a + b) % 2 == 0 else -1
            v = (Qc(a, b) + Qc(a, -b) + Qc(-a, b) + Qc(-a, -b)).real * sg / 2
            if a == 0 and b == 0:
                v = Qc(0, 0).real
            elif a == 0 or b == 0:
                v = v / r2
            Qe[a, b] = v.mid(); Qe[b, a] = v.mid()
            if a >= 1 and b >= 1:
                v = (Qc(a, b) - Qc(a, -b) - Qc(-a, b) + Qc(-a, -b)).real * sg / 2
                Qo[a - 1, b - 1] = v.mid(); Qo[b - 1, a - 1] = v.mid()
    return Qe, Qo


# ---------------- prolates and the E-map ----------------
def prolates(lam2, count):
    c = 2 * math.pi * lam2
    kmax = int(2 * c + 200)
    ks = list(range(0, kmax + 1, 2))
    n = len(ks)
    cA = 2 * PI * lam2
    c2 = cA * cA
    diag = [arb(k * (k + 1)) + c2 * (2 * k * k + 2 * k - 1) / arb((2 * k - 1) * (2 * k + 3)) for k in ks]
    off = [c2 * (k + 1) * (k + 2) / (arb(2 * k + 3) * arb((2 * k + 1) * (2 * k + 5)).sqrt()) for k in ks[:-1]]
    # double-precision eigenvalues as shifts, then inverse iteration in Arb (tridiagonal solves)
    D = np.array([float(d.mid()) for d in diag]); O = np.array([float(o.mid()) for o in off])
    wv = np.linalg.eigvalsh(np.diag(D) + np.diag(O, 1) + np.diag(O, -1))[:count]
    vecs, chis = [], []
    for j in range(count):
        shift = arb(repr(float(wv[j])))
        v = [ONE] * n
        for it in range(10):
            # solve (T - shift) x = v, Thomas algorithm
            a = [diag[i] - shift for i in range(n)]
            cp = [ZERO] * n; dp = [ZERO] * n
            cp[0] = off[0] / a[0]; dp[0] = v[0] / a[0]
            for i in range(1, n):
                den = a[i] - off[i - 1] * cp[i - 1]
                cp[i] = off[i] / den if i < n - 1 else ZERO
                dp[i] = (v[i] - off[i - 1] * dp[i - 1]) / den
            x = [ZERO] * n; x[-1] = dp[-1]
            for i in range(n - 2, -1, -1):
                x[i] = dp[i] - cp[i] * x[i + 1]
            nrm = sum((xi * xi for xi in x), ZERO).sqrt()
            v = [(xi / nrm).mid() for xi in x]
        if float(v[0].mid()) < 0:
            v = [-xi for xi in v]
        vecs.append(v)
    return ks, vecs


def acl_vectors(lam2, N, verbose=True):
    t0 = time.time()
    lam = arb(lam2).sqrt(); L = arb(lam2).log()
    count = 2 * lam2 - 1
    ks, B = prolates(lam2, count)
    nk = len(ks)
    norm = [arb((2 * k + 1)).sqrt() / arb(2).sqrt() for k in ks]      # Pbar_k = sqrt((2k+1)/2) P_k
    I = [lam * arb(2).sqrt() * B[j][0] for j in range(count)]
    # value at 0: P_k(0)
    vals0 = [ONE, ZERO]
    for k in range(1, ks[-1]):
        vals0.append(-k * vals0[k - 1] / (k + 1))
    h0 = [sum((B[j][i] * norm[i] * vals0[ks[i]] for i in range(nk)), ZERO) for j in range(count)]
    kappa = [abs(I[j] / h0[j]) for j in range(count)]
    if verbose:
        print(f"  prolates: {count} states, {nk} Legendre terms ({time.time() - t0:.0f} s)"); sys.stdout.flush()
    # t-quadrature split at the jumps t = log(lam/m)
    jumps = [L - arb(m).log() for m in range(2, lam2)]                 # s = L - log m: E-map jumps
    wN = 2 * math.pi * N / float(L.mid())
    tn = quad_nodes(L, jumps, wN, grade=False)
    Bmat = arb_mat(nk, count)
    for j in range(count):
        for i in range(nk):
            Bmat[i, j] = B[j][i] * norm[i]
    E = [[ZERO] * count for _ in tn]
    rows, owner, npairs = [], [], 0
    def flush():
        if not rows: return
        Y = arb_mat(rows) * Bmat
        for r in range(Y.nrows()):
            o = owner[r]
            for j in range(count):
                E[o][j] += Y[r, j]
        rows.clear(); owner.clear()
    for idx, (s, w) in enumerate(tn):
        t = s - L / 2
        u = t.exp()
        mmax = min(int(math.floor(float((lam / u).mid()) + 1e-12)), lam2)
        for m in range(1, mmax + 1):
            x = m * u / lam
            P = [ONE, x]
            for k in range(1, ks[-1]):
                P.append(((2 * k + 1) * x * P[k] - k * P[k - 1]) / (k + 1))
            rows.append([P[k] for k in ks]); owner.append(idx); npairs += 1
            if len(rows) >= 1500:
                flush()
    flush()
    if verbose:
        print(f"  E-map at {len(tn)} nodes, {npairs} (node, m) pairs ({time.time() - t0:.0f} s)"); sys.stdout.flush()
    # projections on xi_j
    sqrt2L = (2 / L).sqrt(); invsqL = 1 / L.sqrt()
    Xc = arb_mat(N + 1, len(tn)); Xs = arb_mat(N, len(tn)); Em = arb_mat(len(tn), count)
    w1 = 2 * PI / L
    for idx, (s, w) in enumerate(tn):
        t = s - L / 2
        su = (t / 2).exp()
        for j in range(count):
            Em[idx, j] = (E[idx][j] * su).mid()
        s1, c1 = (w1 * t).sin_cos()
        sm, cm = ZERO, ONE
        for m in range(N + 1):
            Xc[m, idx] = (w * (invsqL if m == 0 else sqrt2L) * cm).mid()
            if m >= 1:
                Xs[m - 1, idx] = (w * sqrt2L * sm).mid()
            sm, cm = sm * c1 + cm * s1, cm * c1 - sm * s1
    Pc = Xc * Em; Ps = Xs * Em
    Ve = arb_mat(N + 1, lam2); Vo = arb_mat(N, lam2 - 1)
    for k in range(1, lam2):
        for m in range(N + 1):
            Ve[m, k - 1] = (I[0] * Pc[m, 2 * k] - I[2 * k] * Pc[m, 0]).mid()
    for k in range(1, lam2 - 1):
        for m in range(N):
            Vo[m, k - 1] = (I[1] * Ps[m, 2 * k + 1] - I[2 * k + 1] * Ps[m, 1]).mid()
    # pole vectors (exact): rho_ev = 2 cosh(t/2), rho_odd = 2 sinh(t/2)
    a = HALF; sh = (a * L / 2).sinh()
    for m in range(N + 1):
        om = 2 * PI * m / L; sg = 1 if m % 2 == 0 else -1
        val = 2 * 2 * a * sh * sg / (a * a + om * om)
        Ve[m, lam2 - 1] = (val * (invsqL if m == 0 else sqrt2L)).mid()
        if m >= 1:
            Vo[m - 1, lam2 - 2] = (-2 * 2 * om * sh * sg / (a * a + om * om) * sqrt2L).mid()
    d_ev = [1 - kappa[2 * k] for k in range(1, lam2)]
    d_odd = [1 - kappa[2 * k + 1] for k in range(1, lam2 - 1)]
    if verbose:
        print(f"  projections done ({time.time() - t0:.0f} s)"); sys.stdout.flush()
    return dict(ev=Ve, odd=Vo, d_ev=d_ev, d_odd=d_odd)


# ---------------- linear algebra ----------------
def gram_schmidt(cols, n):
    """Orthonormalize the given column vectors (lists of arb) twice; returns list of vectors."""
    out = []
    for v in cols:
        for _ in range(2):
            for q in out:
                d = sum((vi * qi for vi, qi in zip(v, q)), ZERO)
                v = [vi - d * qi for vi, qi in zip(v, q)]
        nrm = sum((vi * vi for vi in v), ZERO).sqrt()
        out.append([(vi / nrm).mid() for vi in v])
    return out


def bases(V):
    n, k = V.nrows(), V.ncols()
    cols = [[V[i, j] for i in range(n)] for j in range(k)]
    Vq = gram_schmidt(cols, n)
    # complement: project the unit vectors, keep the n - k largest residuals greedily
    W = []
    basis = list(Vq)
    for i in range(n):
        e = [ONE if r == i else ZERO for r in range(n)]
        for _ in range(2):
            for q in basis:
                d = q[i] if _ == 0 and False else sum((ei * qi for ei, qi in zip(e, q)), ZERO)
                e = [ei - d * qi for ei, qi in zip(e, q)]
        nrm = sum((ei * ei for ei in e), ZERO).sqrt()
        if float(nrm.mid()) > 1e-8:
            w = [(ei / nrm).mid() for ei in e]
            basis.append(w); W.append(w)
        if len(W) == n - k:
            break
    Vm = arb_mat(n, k); Wm = arb_mat(n, len(W))
    for j, q in enumerate(Vq):
        for i in range(n): Vm[i, j] = q[i]
    for j, q in enumerate(W):
        for i in range(n): Wm[i, j] = q[i]
    return Vm, Wm


def schur_parts(Q, Vm, Wm):
    A = Vm.transpose() * Q * Vm
    Bm = Vm.transpose() * Q * Wm
    D = Wm.transpose() * Q * Wm
    X = arb_mat.solve(D, Bm.transpose())          # D^{-1} B^T
    S = A - Bm * X
    return A, Bm, D, X, S


def to_mp(M):
    return mpmath.matrix([[mpmath.mpf(M[i, j].mid().str(DIG, radius=False)) for j in range(M.ncols())] for i in range(M.nrows())])


def sym_eig(M):
    Mm = to_mp(M)
    Mm = (Mm + Mm.T) / 2
    E, Q = mpmath.eigsy(Mm)
    return E, Q
