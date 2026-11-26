"""Is the one state that falls when the primes are zoomed in the pole direction? (notes/questions.md)
Pole vectors as in Andrews-Connes-La Gatta: 2 cosh(t/2) (even), 2 sinh(t/2) (odd) on the window, in the Fourier basis.
For each move: negative eigenvalues of the form, the falling state's overlap with the pole vector, and the negative
eigenvalues after restricting to functions orthogonal to the pole vector.
Usage: python3 zoom_pole.py N c p delta...    (p = 0: zoom p -> p^(1 + delta); else that prime moved to p (1 + delta))"""
import sys, math
import numpy as np
from flint import arb
import hp_block as hb
N = int(sys.argv[1]); c = float(sys.argv[2]); p0 = int(sys.argv[3]); F = hb.Form(c); L = F.L; Lf = float(L.mid()); primes = sorted({pp for q, pp in F.pp}); ALL = tuple(q for q, pp in F.pp)
m = np.arange(N + 1); om = 2 * math.pi * m / Lf; sg = (-1.0) ** m; sh = math.sinh(Lf / 4); nrm = np.where(m == 0, 1 / math.sqrt(Lf), math.sqrt(2 / Lf))
pole = [2 * sh * sg / (0.25 + om ** 2) * nrm, (-4 * om * sh * sg / (0.25 + om ** 2) * nrm)[1:]]; pole = [v / np.linalg.norm(v) for v in pole]
def form(delta):
    e = arb(str(delta)); xs = [arb(p) ** (1 + e) if p0 == 0 else (arb(p) * (1 + e) if p == p0 else arb(p)) for p in primes]
    terms = [(x.log(), x ** k) for x in xs for k in range(1, 40) if float((x ** k).mid()) < c]
    Hf = lambda U: F.H(U, exclude=ALL) + sum((lg * F.ramp(U, q) for lg, q in terms), hb.ZERO)
    nodes = hb.quad_nodes(L, [q.log() for lg, q in terms], 2 * math.pi * N / Lf)
    S_, C2_ = hb.transforms(nodes, [Hf(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, Hf(hb.ZERO), Hf(L), L, N); out = []
    for Q, n in ((Qe, N + 1), (Qo, N)):
        M = np.array([[float(Q[i, j].mid()) for j in range(n)] for i in range(n)]); out.append((M + M.T) / 2)
    return out
for delta in [float(x) for x in sys.argv[4:]]:
    print(f"   delta {delta:+.0e}:")
    for h, name in ((0, "even"), (1, "odd ")):
        M = form(delta)[h] if h == 0 else Ms[1]
        if h == 0: Ms = form(delta); M = Ms[0]
        E, U = np.linalg.eigh(M); neg = E[E < -1e-9]; v = pole[h]; B = np.linalg.svd(np.eye(len(v)) - np.outer(v, v))[0][:, :len(v) - 1]; E2 = np.linalg.eigvalsh(B.T @ M @ B); neg2 = E2[E2 < -1e-9]
        print(f"      {name}: {len(neg)} negative (lowest {E[0]:+.2e}); falling state's overlap with the pole vector {abs(U[:, 0] @ v):.3f}; without the pole direction: {len(neg2)} negative (lowest {E2[0]:+.2e})", flush=True)
