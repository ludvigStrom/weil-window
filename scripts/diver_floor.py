"""How far does the form fall when one prime's terms are left out? (notes/function_field.md)
Function-field toy: without one point of degree d the lowest eigenvalue tends to the minimum of the removed term's symbol,
-2 d r/(1 + r), r = q^(-d/2). The same computation for a prime p (all its powers left out) gives -2 ln p/(sqrt p + 1).
Odd and even halves, double-precision eigenvalues of the ball-arithmetic form.  Usage: python3 diver_floor.py N p c1 c2 ..."""
import sys, math
import numpy as np
from flint import arb
import hp_block as hb
N = int(sys.argv[1]); ps = [int(x) for x in sys.argv[2].split(",")]; print(f"primes {ps} left out (all their powers); long-run floors -2 ln p/(sqrt p + 1): " + ", ".join(f"{-2 * math.log(p) / (math.sqrt(p) + 1):.4f}" for p in ps) + f"; first-power floors -ln p/sqrt p: " + ", ".join(f"{-math.log(p) / math.sqrt(p):.4f}" for p in ps))
for c in [float(x) for x in sys.argv[3:]]:
    F = hb.Form(c); L = F.L; EX = tuple(q for q, pp in F.pp if pp in ps)
    nodes = hb.quad_nodes(L, [arb(q).log() for q, pp in F.pp if q < c], 2 * math.pi * N / float(L.mid()))
    Hv = [F.H(s, exclude=EX) for s, w in nodes]; S_, C2_ = hb.transforms(nodes, Hv, L, N); Qe, Qo = hb.gram(S_, C2_, F.H(hb.ZERO, exclude=EX), F.H(L, exclude=EX), L, N)
    lo = []
    for Q, n in ((Qe, N + 1), (Qo, N)):
        M = np.array([[float(Q[i, j].mid()) for j in range(n)] for i in range(n)]); lo.append(np.linalg.eigvalsh((M + M.T) / 2)[0])
    print(f"   c = {c:6.2f}: lowest eigenvalue, even half {lo[0]:+.4f}, odd half {lo[1]:+.4f}", flush=True)
