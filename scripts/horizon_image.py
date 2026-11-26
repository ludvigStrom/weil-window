"""Inner horizons (notes/picture.md, round 5, test 1).
In the window [-a, a] (x = log u) the copy of the window shrunk by n has its edge at depth w_n = L - log n below the top.
The eigen-equation Q u = eps u contains lambda_n u(x - log n), which switches on at that depth carrying the edge profile of u at the
bottom end: -lambda_n beta / sqrt(ln(1/(kappa s))) for s = w_n - t > 0 (odd states). If the local operator is ln(|D|/2 pi), the state
should carry an image of its own edge there:  u_image(s) = -sign * lambda_n beta / ln(1/(kappa s))^(3/2) for s > 0, nothing for s < 0,
so the slope has a one-sided spike  lambda_n beta * 1.5 / (s ln(1/(kappa s))^(5/2)).
Printed: the slope du/dt on both sides of each inner horizon, the difference, and the predicted spike.
Usage: python3 horizon_image.py N c K"""
import sys, math, time, os
import numpy as np, mpmath
from mpmath import mp
from flint import arb, arb_mat
import hp_block as hb
mp.prec = hb.PREC
N = int(sys.argv[1]); c0 = float(sys.argv[2]); K = int(sys.argv[3]); t0 = time.time(); KAP = 2 * math.pi * math.exp(0.5772156649)
def tomp(x): return mp.mpf(x.mid().str(hb.DIG, radius=False))
def lowest(c, vec=False):
    F = hb.Form(c); L = F.L
    nodes = hb.quad_nodes(L, [arb(q).log() for q, p in F.pp if q < c], 2 * math.pi * N / float(L.mid()))
    Hv = [F.H(s) for s, w in nodes]; S_, C2_ = hb.transforms(nodes, Hv, L, N); Qo = hb.gram(S_, C2_, F.H(hb.ZERO), F.H(L), L, N)[1]
    V = arb_mat(N, K, [arb(math.sin(1.7 * (i + 1) * (j + 1)) / (i + 1)) for i in range(N) for j in range(K)])
    for it in range(8):
        V = Qo.solve(V); cols = hb.gram_schmidt([[V[i, j] for i in range(N)] for j in range(K)], N); V = arb_mat(N, K, [cols[j][i] for i in range(N) for j in range(K)])
    H = V.transpose() * Qo * V; Hm = mpmath.matrix(K, K)
    for i in range(K):
        for j in range(K): Hm[i, j] = (tomp(H[i, j]) + tomp(H[j, i])) / 2
    E, W = mpmath.eigsy(Hm); order = sorted(range(K), key=lambda i: E[i]); U = V * arb_mat(K, K, [arb(mpmath.nstr(W[i, order[j]], hb.DIG)) for i in range(K) for j in range(K)])
    return [E[i] for i in order], (np.array([[float(U[i, j].mid()) for j in range(K)] for i in range(N)]) if vec else None), float(L.mid()), F
dc = 0.004 * c0; E0, Uf, Lf, F = lowest(c0, True); Em = lowest(c0 - dc)[0]; Ep = lowest(c0 + dc)[0]
da = (math.log(c0 + dc) - math.log(c0 - dc)) / 2; beta = [math.sqrt(max(float(-(Ep[j] - Em[j]) / da / 2), 0.0)) for j in range(K)]
om = 2 * math.pi * (np.arange(N) + 1) / Lf
def u(t): return math.sqrt(2 / Lf) * np.sin(om * (Lf / 2 - t)) @ Uf
def du(t): return -math.sqrt(2 / Lf) * (om * np.cos(om * (Lf / 2 - t))) @ Uf
sgn = np.sign(u(0.004))          # sign of each state just below the top edge
print(f"c = {c0}, N = {N}, L = {Lf:.4f}  ({time.time() - t0:.0f} s)"); print("  j | eps | beta"); [print(f"  {j + 1} | {float(E0[j]):.3e} | {beta[j]:.3e}") for j in range(K)]
np.savez(f"horizon_image_{c0}_{N}.npz", U=Uf, L=Lf, beta=np.array(beta), eps=np.array([float(e) for e in E0]))
# detector: D(s) = u(w - s) - 2 u(w) + u(w + s); smooth background is ~ s^2, so E(s) = [4 D(s/2) - D(s)] / 3 removes it and keeps
# the one-sided image g(s) = -sign lambda beta / ln(1/(kappa s))^1.5:  E = [4 g(s/2) - g(s)] / 3.
def g(s, lam): return lam / math.log(1 / (KAP * s)) ** 1.5
for q, p in F.pp:
    if q >= c0: continue
    w = Lf - math.log(q); lam = math.log(p) / math.sqrt(q)
    if w > Lf / 2 or w < 0.03: continue
    print(f"  inner horizon of n = {q} at depth {w:.4f} (lambda = {lam:.3f}).   s | E(s) / beta for each state, sign-adjusted so that the predicted image is negative | predicted")
    for s in (0.012, 0.016, 0.024, 0.032):
        D = lambda h: (u(w - h) - 2 * u(w) + u(w + h)) * sgn
        E = (4 * D(s / 2) - D(s)) / 3; pred = -(4 * g(s / 2, lam) - g(s, lam)) / 3
        print(f"     {s:.3f} | " + "  ".join(f"{E[j] / beta[j]:9.3f}" for j in range(K)) + f" | {pred:.3f}")
# composite n (no term of its own in the form): is there a second-order image?
for q in [int(x) for x in os.environ.get("EXTRA", "").split(",") if x]:
    w = Lf - math.log(q); print(f"  composite n = {q} (no term in the form), depth {w:.4f}:   s | E(s) / beta")
    for s in (0.012, 0.016, 0.024, 0.032):
        D = lambda h: (u(w - h) - 2 * u(w) + u(w + h)) * sgn; E = (4 * D(s / 2) - D(s)) / 3
        print(f"     {s:.3f} | " + "  ".join(f"{E[j] / beta[j]:9.3f}" for j in range(K)))
# control: the same detector at depths where there is no inner horizon
for w in (0.2, 0.45):
    print(f"  control, depth {w} (no inner horizon)")
    for s in (0.012, 0.024):
        D = lambda h: (u(w - h) - 2 * u(w) + u(w + h)) * sgn; E = (4 * D(s / 2) - D(s)) / 3
        print(f"     {s:.3f} | " + "  ".join(f"{E[j] / beta[j]:9.3f}" for j in range(K)))
