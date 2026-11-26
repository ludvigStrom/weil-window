"""Question 9: what happens to the form when one prime is moved? (notes/questions.md)
The prime p is replaced by a 'prime' at the real position x = p (1 + delta), with weight ln x and powers x^k
(a Beurling-style prime); everything else is unchanged. Prints the lowest eigenvalue of the even and odd halves.
Double-precision eigenvalues of the ball-arithmetic form, so shifts below about 1e-12 are not resolved.
Usage: python3 prime_move.py N c p delta1 delta2 ...   (appends to prime_move.json)"""
import sys, math, json, os
import numpy as np
from flint import arb
import hp_block as hb
N = int(sys.argv[1]); c = float(sys.argv[2]); p = int(sys.argv[3]); out = []
F = hb.Form(c); L = F.L; EX = tuple(q for q, pp in F.pp if pp == p); others = [arb(q).log() for q, pp in F.pp if pp != p]
for delta in [float(x) for x in sys.argv[4:]]:
    x = arb(p) * (1 + arb(str(delta))); pw = [x ** k for k in range(1, 40) if float((x ** k).mid()) < c]
    H = lambda U: F.H(U, exclude=EX) + sum((x.log() * F.ramp(U, q) for q in pw), hb.ZERO)
    nodes = hb.quad_nodes(L, others + [q.log() for q in pw], 2 * math.pi * N / float(L.mid()))
    S_, C2_ = hb.transforms(nodes, [H(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, H(hb.ZERO), H(L), L, N); lo = []
    for Q, n in ((Qe, N + 1), (Qo, N)):
        M = np.array([[float(Q[i, j].mid()) for j in range(n)] for i in range(n)]); lo.append(float(np.linalg.eigvalsh((M + M.T) / 2)[0]))
    out.append(dict(c=c, p=p, N=N, delta=delta, even=lo[0], odd=lo[1]))
    print(f"   prime {p} at {p * (1 + delta):.5f} (delta {delta:+.0e}): lowest eigenvalue, even half {lo[0]:+.3e}, odd half {lo[1]:+.3e}   per unit delta: {lo[0] / delta if delta else 0:+.4f}, {lo[1] / delta if delta else 0:+.4f}", flush=True)
R = json.load(open("prime_move.json")) if os.path.exists("prime_move.json") else []; json.dump(R + out, open("prime_move.json", "w"))
