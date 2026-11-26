"""The strength of each prime by first-order perturbation theory (notes/questions.md, 'The map'), for larger windows.
One high-precision build of the even form Q, its inverse, the vacuum v0 by inverse iteration, then for each prime the
exact derivative P = dQ/d(delta) from the shift operators A_a in closed form:
   dQ/d(delta) = - sum_k p^(-k/2) [ (1 - (k/2) ln p) A_(k ln p) + k ln p A'_(k ln p) ],
the response dv = -(Q^-1 off v0) P v0, and the shift of the first tick  -dv^(t0) / v0^'(t0).
Agrees with sens_map.py where both run.   Usage: python3 sens_fast.py N c p1,p2,...    env DPS (default 170)"""
import sys, math, os, json, time
import mpmath
from mpmath import mp
from flint import arb
import hp_block as hb
mp.dps = int(os.environ.get("DPS", 170)); N = int(sys.argv[1]); c = float(sys.argv[2]); T0 = time.time()
F = hb.Form(c); L = F.L; tomp = lambda x: mp.mpf(x.mid().str(hb.DIG, radius=False)); Lm = tomp(L); n1 = N + 1
om = [2 * mp.pi * m / Lm for m in range(2 * N + 2)]; nrm = [1 / mp.sqrt(Lm)] + [mp.sqrt(2 / Lm)] * N
todo = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else sorted({pp for q, pp in F.pp})
nodes = hb.quad_nodes(L, [arb(q).log() for q, pp in F.pp], 2 * math.pi * N / float(L.mid()))
S_, C2_ = hb.transforms(nodes, [F.H(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, F.H(hb.ZERO), F.H(L), L, N)
Q = mp.matrix(n1, n1)
for i in range(n1):
    for j in range(n1): Q[i, j] = (tomp(Qe[i, j]) + tomp(Qe[j, i])) / 2
print(f"c = {c}, N = {N}: form built ({time.time() - T0:.0f} s)", flush=True); Qi = mp.inverse(Q); v = mp.matrix([mp.mpf(1) / (i + 1) for i in range(n1)])
for it in range(6): v = Qi * v; v = v / mp.norm(v)
lam0 = (v.T * Q * v)[0]; w2 = Qi * v; print(f"   vacuum eigenvalue {mp.nstr(lam0, 4)} (check: 1 / |Q^-1 v| = {mp.nstr(1 / mp.norm(w2), 4)}); inverse done ({time.time() - T0:.0f} s)", flush=True)
def hat(vec):
    def f(s):
        return sum(vec[m] * nrm[m] * (mp.sincpi((s - om[m]) * Lm / (2 * mp.pi)) + mp.sincpi((s + om[m]) * Lm / (2 * mp.pi))) for m in range(n1))
    return f
f0 = hat(v); t0 = mp.findroot(f0, mp.mpf("14.1347"), solver="secant"); d0 = mp.diff(f0, t0); print(f"   first tick {mp.nstr(t0, 15)}; minus the first zero: {mp.nstr(t0 - mp.im(mp.zetazero(1)), 3)}", flush=True)
def Aops(a):
    """A_a and dA_a/da on the even basis, closed form."""
    lo, hi = -Lm / 2, Lm / 2 - a; G = mp.matrix(n1, n1); D = mp.matrix(n1, n1)
    sl = [mp.sin(o * lo) for o in om]; cl = [mp.cos(o * lo) for o in om]; sh = [mp.sin(o * hi) for o in om]; ch = [mp.cos(o * hi) for o in om]; sa = [mp.sin(om[n] * a) for n in range(n1)]; ca = [mp.cos(om[n] * a) for n in range(n1)]; cL = [mp.cos(om[n] * Lm / 2) for n in range(n1)]
    def tab(k):                                    # sin, cos of om_k * x at x = hi, lo for signed index k
        s = 1 if k >= 0 else -1; k = abs(k); return s * sh[k], ch[k], s * sl[k], cl[k]
    for m in range(n1):
        for n in range(n1):
            tot_c = tot_s = mp.mpf(0)
            for k, sgn in ((m - n, -1), (m + n, 1)):                # integrand cos / sin (om_k u + sgn om_n a)
                shk, chk, slk, clk = tab(k); sp, cp = sgn * sa[n], ca[n]
                if k == 0: Ic = (hi - lo) * cp; Is = (hi - lo) * sp
                else:
                    kap = om[abs(k)] * (1 if k > 0 else -1)
                    Ic = ((shk * cp + chk * sp) - (slk * cp + clk * sp)) / kap                     # int cos(kap u + phi)
                    Is = -((chk * cp - shk * sp) - (clk * cp - slk * sp)) / kap                    # int sin(kap u + phi)
                tot_c += Ic; tot_s += sgn * Is
            G[m, n] = nrm[m] * nrm[n] * tot_c / 2; D[m, n] = nrm[m] * nrm[n] * (-om[n] * tot_s / 2 - ch[m] * cL[n])
    return G + G.T, D + D.T
out = {}
for p in todo:
    Pv = mp.matrix(n1, 1); k = 1
    while p ** k < c:
        a = k * mp.log(p); A, dA = Aops(a); Pv += -mp.mpf(p) ** (-mp.mpf(k) / 2) * ((1 - mp.mpf(k) / 2 * mp.log(p)) * (A * v) + k * mp.log(p) * (dA * v)); k += 1
    r = Pv - v * (v.T * Pv)[0]; dv = -(Qi * r); dv = dv - v * (v.T * dv)[0]; S = -hat(dv)(t0) / d0; out[p] = float(mp.log10(abs(S)))
    print(f"   prime {p}: log10 |d(first tick)/d(delta)| = {out[p]:+.2f}   ({time.time() - T0:.0f} s)", flush=True)
fn = f"sens_fast_{c:g}.json"; old = json.load(open(fn)) if os.path.exists(fn) else {}; old.update({str(k): x for k, x in out.items()}); json.dump(old, open(fn, "w"))
