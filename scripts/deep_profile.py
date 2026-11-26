"""Shape of the near-null odd states inside the window (notes/margin.md, part 10).
The K lowest odd eigenvectors by block inverse iteration; values u_j(l - t) at depth t from the right edge, saved to JSON.
Usage: python3 deep_profile.py N c K"""
import sys, math, json, time
import numpy as np, mpmath
from mpmath import mp
from flint import arb, arb_mat
import hp_block as hb
mp.prec = hb.PREC
N = int(sys.argv[1]); c0 = float(sys.argv[2]); K = int(sys.argv[3]); t0 = time.time()
F = hb.Form(c0); L = F.L
nodes = hb.quad_nodes(L, [arb(q).log() for q, p in F.pp if q < c0], 2 * math.pi * N / float(L.mid()))
Hv = [F.H(s) for s, w in nodes]; S_, C2_ = hb.transforms(nodes, Hv, L, N); Qo = hb.gram(S_, C2_, F.H(hb.ZERO), F.H(L), L, N)[1]
def tomp(x): return mp.mpf(x.mid().str(hb.DIG, radius=False))
V = arb_mat(N, K, [arb(math.sin(1.7 * (i + 1) * (j + 1)) / (i + 1)) for i in range(N) for j in range(K)])
for it in range(9):
    V = Qo.solve(V); cols = hb.gram_schmidt([[V[i, j] for i in range(N)] for j in range(K)], N); V = arb_mat(N, K, [cols[j][i] for i in range(N) for j in range(K)])
H = V.transpose() * Qo * V; Hm = mpmath.matrix(K, K)
for i in range(K):
    for j in range(K): Hm[i, j] = (tomp(H[i, j]) + tomp(H[j, i])) / 2
E, W = mpmath.eigsy(Hm); order = sorted(range(K), key=lambda i: E[i])
U = V * arb_mat(K, K, [arb(mpmath.nstr(W[i, order[j]], hb.DIG)) for i in range(K) for j in range(K)])
Uf = np.array([[float(U[i, j].mid()) for j in range(K)] for i in range(N)]); eps = [float(E[i]) for i in order]; Lf = float(L.mid()); om = 2 * math.pi * (np.arange(N) + 1) / Lf
ts = np.concatenate([np.arange(1, 40) * 0.0005, np.arange(1, 200) * Lf / 400 ]); ts = np.unique(np.round(ts[ts < Lf / 2], 6))
vals = np.array([math.sqrt(2 / Lf) * np.sin(om * (Lf / 2 - t)) @ Uf for t in ts])
# spectral content: mean frequency and the frequency below which 90% of the weight lies
P = Uf ** 2; mean = (P * om[:, None]).sum(0); cum = np.cumsum(P, 0); f90 = [float(om[np.searchsorted(cum[:, j], 0.9)]) for j in range(K)]
json.dump(dict(N=N, c=c0, L=Lf, eps=eps, t=ts.tolist(), u=vals.tolist(), mean_freq=mean.tolist(), f90=f90), open(f"deep_profile_{c0}_{N}.json", "w"))
np.savez(f"deep_profile_{c0}_{N}.npz", U=Uf, L=Lf, eps=np.array(eps))
print(f"c = {c0}, N = {N}, K = {K}, L/2 = {Lf / 2:.4f}  ({time.time() - t0:.0f} s)")
print("  j | eps | mean frequency | 90% of weight below | sign changes on the half window")
for j in range(K): print(f"  {j + 1} | {eps[j]:.3e} | {mean[j]:7.1f} | {f90[j]:7.1f} | {int((np.diff(np.sign(vals[ts > 0.004, j])) != 0).sum())}")
