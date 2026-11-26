"""Testing the 'jammed packing' picture (notes/questions.md). Two signatures:
(a) one contact removed -> exactly one mode goes soft: count the negative eigenvalues when prime powers are left out;
(b) soft modes of a jammed solid have a flat density in frequency (number below stiffness k grows like sqrt k):
    how are the tiny levels distributed?   Usage: python3 jam_test.py N c"""
import sys, math, json
import numpy as np
from flint import arb
import hp_block as hb
N = int(sys.argv[1]); c = float(sys.argv[2]); F = hb.Form(c); L = F.L; pp = F.pp
def neg(EX):
    nodes = hb.quad_nodes(L, [arb(q).log() for q, p in pp], 2 * math.pi * N / float(L.mid()))
    S_, C2_ = hb.transforms(nodes, [F.H(s, exclude=EX) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, F.H(hb.ZERO, exclude=EX), F.H(L, exclude=EX), L, N); out = []
    for Q, n in ((Qe, N + 1), (Qo, N)):
        M = np.array([[float(Q[i, j].mid()) for j in range(n)] for i in range(n)]); E = np.linalg.eigvalsh((M + M.T) / 2); out.append(E[E < -1e-9])
    return out
print(f"c = {c}: prime powers in the window: {[q for q, p in pp]}")
cases = [((q,), f"only {q} left out") for q, p in pp] + [(tuple(q for q, p in pp if p == r), f"prime {r} with all its powers") for r in sorted({p for q, p in pp}) if sum(1 for q, p in pp if p == r) > 1]
cases += [((5, 7), "5 and 7"), ((7, 11, 13), "7, 11 and 13"), (tuple(q for q, p in pp), "everything")]
for EX, name in cases:
    e, o = neg(EX); print(f"   {name:32s}: {len(EX)} removed -> {len(e) + len(o)} negative (even {len(e)}, odd {len(o)}); lowest {min(list(e) + list(o) + [0]):+.3f}", flush=True)
