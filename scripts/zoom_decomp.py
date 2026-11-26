"""Where does the cost of the zoom p -> p^(1+eps) come from? (notes/questions.md)
Exact bookkeeping: stretching the test functions turns the window c into c' = c^(1/(1+eps)), and the true form on c'
is positive. So   Q_zoomed(c) = Q_true(c') + eps * D,   and the whole cost is in D. D splits into
  D_smooth : the same difference with every prime removed (archimedean and pole terms, which do not stretch)
  D_prime  : the rest; the prime terms stretch except for their weights, (1 + eps) n^(-eps/2) instead of 1.
Prints the lowest eigenvalue of Q_true(c') + eps * X for X = D, D_smooth, D_prime (per unit eps), even and odd halves.
Usage: python3 zoom_decomp.py N c [eps]"""
import sys, math
import numpy as np
from flint import arb
import hp_block as hb
N = int(sys.argv[1]); c = float(sys.argv[2]); eps = float(sys.argv[3]) if len(sys.argv) > 3 else 1e-3
def mats(c, stretch, with_primes):
    F = hb.Form(c if stretch else c ** (1 / (1 + eps))); L = F.L; primes = sorted({pp for q, pp in F.pp}); ALL = tuple(q for q, pp in F.pp); terms = []
    if with_primes:
        e = arb(str(eps)) if stretch else arb(0)
        terms = [(x.log(), x ** k) for x in (arb(p) ** (1 + e) for p in primes) for k in range(1, 40) if float((x ** k).mid()) < float(F.lam2)]
    Hf = lambda U: F.H(U, exclude=ALL) + sum((lg * F.ramp(U, q) for lg, q in terms), hb.ZERO)
    nodes = hb.quad_nodes(L, [q.log() for lg, q in terms], 2 * math.pi * N / float(L.mid()))
    S_, C2_ = hb.transforms(nodes, [Hf(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, Hf(hb.ZERO), Hf(L), L, N); out = []
    for Q, n in ((Qe, N + 1), (Qo, N)):
        M = np.array([[float(Q[i, j].mid()) for j in range(n)] for i in range(n)]); out.append((M + M.T) / 2)
    return out
Zp, Tp, Zn, Tn = mats(c, True, True), mats(c, False, True), mats(c, True, False), mats(c, False, False)
print(f"c = {c}, eps = {eps}: true form on the shrunk window c' = {c ** (1 / (1 + eps)):.4f}")
for h, name in ((0, "even"), (1, "odd")):
    D = (Zp[h] - Tp[h]) / eps; Ds = (Zn[h] - Tn[h]) / eps; Dp = D - Ds; lo = lambda X: np.linalg.eigvalsh(Tp[h] + eps * X)[0] / eps
    E, U = np.linalg.eigh(Tp[h] + eps * D); v = U[:, 0]
    print(f"   {name}: lowest eigenvalue per unit eps: full {lo(D):+.4f} | smooth terms only {lo(Ds):+.4f} | prime weights only {lo(Dp):+.4f}")
    print(f"        on the falling state: true form {v @ Tp[h] @ v / eps:+.4f}, smooth {v @ Ds @ v:+.4f}, prime weights {v @ Dp @ v:+.4f}   (sum = full)")
