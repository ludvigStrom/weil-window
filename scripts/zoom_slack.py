"""Does the smallest eigenvalue measure how much zoom-out the primes tolerate? (notes/questions.md)
Zoom p -> p^(1 + eps). For each half: eigenvalues at eps = 0, the common sinking rate r (from a zoom well above the
smallest eigenvalue), the predicted critical zoom eps* = lowest eigenvalue / its own rate, and the measured eps*
(bisection on the sign of the lowest eigenvalue). Double-precision eigenvalues, so use windows whose smallest
eigenvalue is above 1e-11.   Usage: python3 zoom_slack.py N c"""
import sys, math
import numpy as np
from flint import arb
import hp_block as hb
N = int(sys.argv[1]); c = float(sys.argv[2]); F = hb.Form(c); L = F.L; Lf = float(L.mid()); primes = sorted({pp for q, pp in F.pp}); ALL = tuple(q for q, pp in F.pp)
def form(eps):
    e = arb(f"{float(eps):.17e}"); terms = [(x.log(), x ** k) for x in (arb(p) ** (1 + e) for p in primes) for k in range(1, 40) if float((x ** k).mid()) < c]
    Hf = lambda U: F.H(U, exclude=ALL) + sum((lg * F.ramp(U, q) for lg, q in terms), hb.ZERO)
    nodes = hb.quad_nodes(L, [q.log() for lg, q in terms], 2 * math.pi * N / Lf)
    S_, C2_ = hb.transforms(nodes, [Hf(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, Hf(hb.ZERO), Hf(L), L, N); out = []
    for Q, n in ((Qe, N + 1), (Qo, N)):
        M = np.array([[float(Q[i, j].mid()) for j in range(n)] for i in range(n)]); out.append((M + M.T) / 2)
    return out
M0 = form(0.0); print(f"c = {c}, primes {primes}")
for h, name in ((0, "even"), (1, "odd")):
    E, U = np.linalg.eigh(M0[h]); print(f"   {name}: lowest eigenvalues " + ", ".join(f"{x:.3e}" for x in E[:5]))
    if E[0] < 1e-11: print("      too small for double precision"); continue
    v = U[:, 0]; big = 1e-3; r = -(v @ (form(big)[h] - M0[h]) @ v) / big; pred = E[0] / r     # first-order rate of the lowest state itself
    lo, hi = 0.0, 4 * pred
    while np.linalg.eigvalsh(form(hi)[h])[0] > 0: hi *= 2
    for it in range(18):
        mid = (lo + hi) / 2
        if np.linalg.eigvalsh(form(mid)[h])[0] > 0: lo = mid
        else: hi = mid
    print(f"      sinking rate of the lowest state {r:.4f} per unit zoom; predicted critical zoom-out {pred:.4e}; measured {0.5 * (lo + hi):.4e}; ratio {0.5 * (lo + hi) / pred:.4f}", flush=True)
