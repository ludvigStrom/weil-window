"""Bombieri's test on a moved prime (notes/questions.md): how many negative eigenvalues appear, and how negative?
(Bombieri 2000: with finitely many zeros off the line, the number of negative eigenvalues is half their number.)
Prime p moved to p (1 + delta); or p = 0: all primes zoomed, p -> p^(1 + delta).  Usage: python3 neg_count.py N c p delta..."""
import sys, math
import numpy as np
from flint import arb
import hp_block as hb
N = int(sys.argv[1]); c = float(sys.argv[2]); p0 = int(sys.argv[3]); F = hb.Form(c); L = F.L; primes = sorted({pp for q, pp in F.pp}); ALL = tuple(q for q, pp in F.pp)
for delta in [float(x) for x in sys.argv[4:]]:
    e = arb(str(delta)); xs = [arb(p) ** (1 + e) if p0 == 0 else (arb(p) * (1 + e) if p == p0 else arb(p)) for p in primes]
    terms = [(x.log(), x ** k) for x in xs for k in range(1, 40) if float((x ** k).mid()) < c]
    Hf = lambda U: F.H(U, exclude=ALL) + sum((lg * F.ramp(U, q) for lg, q in terms), hb.ZERO)
    nodes = hb.quad_nodes(L, [q.log() for lg, q in terms], 2 * math.pi * N / float(L.mid()))
    S_, C2_ = hb.transforms(nodes, [Hf(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, Hf(hb.ZERO), Hf(L), L, N); out = []
    for Q, n in ((Qe, N + 1), (Qo, N)):
        M = np.array([[float(Q[i, j].mid()) for j in range(n)] for i in range(n)]); E = np.linalg.eigvalsh((M + M.T) / 2); neg = E[E < -1e-9]
        out.append(f"{len(neg)} negative ({', '.join(f'{x:.1e}' for x in neg[:6])}{', ...' if len(neg) > 6 else ''})")
    print(f"   delta {delta:+.0e}: even {out[0]} | odd {out[1]}", flush=True)
