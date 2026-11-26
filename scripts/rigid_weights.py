"""Number-field rigidity with the weights free as well as the positions (notes/questions.md).
Parameters: for each prime p <= c, its position (p -> p(1 + d_p), from the cache of rigid.py) and its weight
(ln p -> (1 + e_p) ln p on all its powers; the form is linear in e_p, so the derivative is exact).
At fixed positions the set of weights with a positive form is convex; is the true point alone in it?
Usage: python3 rigid_weights.py N c    (needs rigid_P_<c>_<N>.npz; env TOL default 1e-2)"""
import sys, math, os
import numpy as np
from flint import arb
import hp_block as hb, soft
N = int(sys.argv[1]); c = float(sys.argv[2]); TOL = float(os.environ.get("TOL", 1e-2)); Z = dict(np.load(f"rigid_P_{c:g}_{N}.npz"))
F = hb.Form(c); L = F.L; primes = sorted({pp for q, pp in F.pp}); cache = f"rigid_W_{c:g}_{N}.npz"; W = dict(np.load(cache)) if os.path.exists(cache) else {}
def form(EX):
    nodes = hb.quad_nodes(L, [arb(q).log() for q, pp in F.pp], 2 * math.pi * N / float(L.mid()))
    S_, C2_ = hb.transforms(nodes, [F.H(s, exclude=EX) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, F.H(hb.ZERO, exclude=EX), F.H(L, exclude=EX), L, N); out = []
    for Q, n in ((Qe, N + 1), (Qo, N)):
        M = np.array([[float(Q[i, j].mid()) for j in range(n)] for i in range(n)]); out.append((M + M.T) / 2)
    return out
if "Te" not in W: W["Te"], W["To"] = form(()); 
for p in primes:
    if f"We{p}" not in W:
        A = form(tuple(q for q, pp in F.pp if pp == p)); W[f"We{p}"], W[f"Wo{p}"] = W["Te"] - A[0], W["To"] - A[1]; np.savez(cache, **W)   # d/de_p = minus the prime's term
V = []
for M in (Z["Q0e"], Z["Q0o"]):
    E, U = np.linalg.eigh(M); V.append(U[:, E < TOL])
Pw = [[V[0].T @ W[f"We{p}"] @ V[0], V[1].T @ W[f"Wo{p}"] @ V[1]] for p in primes]; Pp = [[V[0].T @ Z[f"De{p}"] @ V[0], V[1].T @ Z[f"Do{p}"] @ V[1]] for p in primes]
rng = np.random.default_rng(3); fmt = lambda d: ", ".join(f"{x:+.2f}" for x in d)
print(f"c = {c}, primes {primes}; near-null space {V[0].shape[1]} + {V[1].shape[1]} states")
print("   one weight alone (per unit relative change, up / down): " + ", ".join(f"{p}: {soft.lam(Pw, np.eye(len(primes))[i]):+.2f}/{soft.lam(Pw, -np.eye(len(primes))[i]):+.2f}" for i, p in enumerate(primes)))
for name, P in (("weights only", Pw), ("positions only", Pp), ("weights and positions", Pw + Pp)):
    s, Fl = soft.flat(P); sp = soft.spectrum(P, rng, starts=8, iters=4000)
    print(f"   {name}: {len(Fl)} flat directions; softness spectrum: " + ", ".join(f"{v:+.3f}" for v, d in sp)); print(f"        softest direction: {fmt(sp[0][1])}")
all_down = -np.ones(len(primes)) / math.sqrt(len(primes)); print(f"   all weights down together (towards 'no primes'): {soft.lam(Pw, all_down):+.3f} per unit;  all up: {soft.lam(Pw, -all_down):+.3f}")
