"""Track the lowest eigenvalues of the even and odd halves of the Weil form as the window c grows (notes/ladder.md).
Block inverse iteration in ball arithmetic. Output: ladder_track.json, entries {c, N, even: [...], odd: [...]}.
Usage: python3 ladder_track.py N K c1 c2 ..."""
import sys, math, json, os, time, fcntl
import mpmath
from mpmath import mp
from flint import arb, arb_mat
import hp_block as hb
mp.prec = hb.PREC
N = int(sys.argv[1]); K = int(sys.argv[2]); cs = [float(x) for x in sys.argv[3:]]
def tomp(x): return mp.mpf(x.mid().str(hb.DIG, radius=False))
def lowest(Q, n):
    V = arb_mat(n, K, [arb(math.sin(1.7 * (i + 1) * (j + 1)) / (i + 1)) for i in range(n) for j in range(K)])
    for it in range(12):
        V = Q.solve(V); cols = hb.gram_schmidt([[V[i, j] for i in range(n)] for j in range(K)], n); V = arb_mat(n, K, [cols[j][i] for i in range(n) for j in range(K)])
    H = V.transpose() * Q * V; Hm = mpmath.matrix(K, K)
    for i in range(K):
        for j in range(K): Hm[i, j] = (tomp(H[i, j]) + tomp(H[j, i])) / 2
    return sorted(float(mpmath.log10(abs(e))) * (1 if e > 0 else -1) if False else float(e) if abs(e) > 1e-300 else float(mpmath.log10(e)) for e in mpmath.eigsy(Hm, eigvals_only=True)), [mpmath.nstr(e, 8) for e in sorted(mpmath.eigsy(Hm, eigvals_only=True))]
for c in cs:
    t0 = time.time(); F = hb.Form(c); L = F.L
    nodes = hb.quad_nodes(L, [arb(q).log() for q, p in F.pp if q < c], 2 * math.pi * N / float(L.mid()))
    EX = tuple(int(x) for x in os.environ.get('EXCL', '').split(',') if x)
    Hv = [F.H(s, exclude=EX) for s, w in nodes]; S_, C2_ = hb.transforms(nodes, Hv, L, N); Qe, Qo = hb.gram(S_, C2_, F.H(hb.ZERO, exclude=EX), F.H(L, exclude=EX), L, N)
    ev = lowest(Qe, N + 1)[1]; od = lowest(Qo, N)[1]
    with open("ladder_track.lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX); f = os.environ.get("OUT", "ladder_track.json"); R = json.load(open(f)) if os.path.exists(f) else []
        R.append(dict(c=c, N=N, even=ev, odd=od)); json.dump(R, open(f, "w"), indent=0)
    print(f"c = {c}: {time.time() - t0:.0f} s; odd {od[:3]} ... {od[-2:]}", flush=True)
