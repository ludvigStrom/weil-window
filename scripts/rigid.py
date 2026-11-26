"""Is the true configuration of primes rigid? (notes/questions.md, rigidity)
Every prime p <= c is moved to p (1 + d_p) (weight ln x, powers x^k). To first order the lowest eigenvalue is the lowest
eigenvalue of sum_p d_p P_p on the near-null space V of the true form (eigenvalues below TOL), P_p = dQ/dd_p restricted to V
(central differences, step H). Rigid to first order iff no direction d makes that matrix positive semidefinite on V.
We maximise lambda_min(sum d_p P_p) over unit vectors d (it is concave, so local search finds the maximum),
then check the best direction on the true form.  Even and odd halves together.
Usage: python3 rigid.py N c     (env H default 1e-4, TOL default 1e-2: the transitional states must be included; writes rigid_<c>.json)"""
import sys, math, json, os
import numpy as np
from scipy.optimize import minimize
from flint import arb
import hp_block as hb
N = int(sys.argv[1]); c = float(sys.argv[2]); H = float(os.environ.get("H", 1e-4)); TOL = float(os.environ.get("TOL", 1e-2))
F = hb.Form(c); L = F.L; primes = sorted({pp for q, pp in F.pp}); ALL = tuple(q for q, pp in F.pp)
def form(d):
    pw = [(arb(p) * (1 + arb(str(d.get(p, 0.0)))), p) for p in primes]; terms = [(x.log(), x ** k) for x, p in pw for k in range(1, 40) if float((x ** k).mid()) < c]
    Hf = lambda U: F.H(U, exclude=ALL) + sum((lg * F.ramp(U, q) for lg, q in terms), hb.ZERO)
    nodes = hb.quad_nodes(L, [q.log() for lg, q in terms], 2 * math.pi * N / float(L.mid()))
    S_, C2_ = hb.transforms(nodes, [Hf(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, Hf(hb.ZERO), Hf(L), L, N); out = []
    for Q, n in ((Qe, N + 1), (Qo, N)):
        M = np.array([[float(Q[i, j].mid()) for j in range(n)] for i in range(n)]); out.append((M + M.T) / 2)
    return out
cache = f"rigid_P_{c:g}_{N}.npz"; Z = dict(np.load(cache)) if os.path.exists(cache) else {}
if "Q0e" not in Z: Z["Q0e"], Z["Q0o"] = form({}); np.savez(cache, **Z)
for p in primes:                                   # resumable: one prime per step, saved as we go
    if f"De{p}" not in Z:
        A, B = form({p: H}), form({p: -H}); Z[f"De{p}"], Z[f"Do{p}"] = (A[0] - B[0]) / (2 * H), (A[1] - B[1]) / (2 * H); np.savez(cache, **Z); print(f"   derivative for prime {p} done", flush=True)
V = []; EV = []
for M in (Z["Q0e"], Z["Q0o"]):
    E, U = np.linalg.eigh(M); V.append(U[:, E < TOL]); EV.append(E[E < TOL])
print(f"c = {c}: primes {primes}; near-null space (eigenvalues below {TOL:g}): even {V[0].shape[1]}, odd {V[1].shape[1]}", flush=True)
P = {p: [V[0].T @ Z[f"De{p}"] @ V[0], V[1].T @ Z[f"Do{p}"] @ V[1]] for p in primes}
def lam(d):
    d = np.asarray(d) / np.linalg.norm(d); return min(np.linalg.eigvalsh(sum(d[i] * P[p][h] for i, p in enumerate(primes)))[0] for h in range(2))
rng = np.random.default_rng(1); best = (-1e9, None)
for t in range(30):
    r = minimize(lambda d: -lam(d), rng.normal(size=len(primes)), method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9, maxiter=4000))
    if -r.fun > best[0]: best = (-r.fun, r.x / np.linalg.norm(r.x))
worst = min(lam(rng.normal(size=len(primes))) for t in range(2000))
print(f"   best direction (unit vector over primes {primes}): " + ", ".join(f"{x:+.3f}" for x in best[1])); print(f"   lowest eigenvalue per unit displacement along it: {best[0]:+.4f}  (a random direction: down to {worst:+.1f}; single primes: " + ", ".join(f"{p}: {lam(np.eye(len(primes))[i]):+.1f}/{lam(-np.eye(len(primes))[i]):+.1f}" for i, p in enumerate(primes)) + ")")
chk = []
for s_ in (1e-3,):
    M = form({p: s_ * best[1][i] for i, p in enumerate(primes)}); lo = [float(np.linalg.eigvalsh(m)[0]) for m in M]; chk.append((s_, lo)); print(f"   true form along the best direction, size {s_:g}: lowest eigenvalue even {lo[0]:+.3e}, odd {lo[1]:+.3e}   (first-order prediction {best[0] * s_:+.3e})", flush=True)
lnp = np.log(primes); lnp /= np.linalg.norm(lnp); print(f"   overlap of the best direction with ln p: {abs(best[1] @ lnp):.4f}")
# the softness spectrum: best direction, then the best one orthogonal to it, and so on
done = []; spec = []
for r_ in range(len(primes)):
    Bc = np.linalg.svd(np.eye(len(primes)) - sum(np.outer(v, v) for v in done))[0][:, :len(primes) - len(done)] if done else np.eye(len(primes)); bb = (-1e9, None)
    if Bc.shape[1] == 1: bb = max((lam(Bc[:, 0]), Bc[:, 0]), (lam(-Bc[:, 0]), -Bc[:, 0]), key=lambda t: t[0])
    else:
        for t in range(25):
            r = minimize(lambda y: -lam(Bc @ y), rng.normal(size=Bc.shape[1]), method="Nelder-Mead", options=dict(xatol=1e-6, fatol=1e-9, maxiter=6000))
            if -r.fun > bb[0]: bb = (-r.fun, Bc @ r.x / np.linalg.norm(r.x))
    done.append(bb[1]); spec.append(bb[0]); print(f"   soft direction {r_ + 1}: {bb[0]:+8.3f} per unit   " + ", ".join(f"{x:+.2f}" for x in bb[1]), flush=True)
# which near-null states hold the pattern: for the zoom direction, the state that goes negative, in the basis of near-null eigenstates
for h, name in ((0, "even"), (1, "odd")):
    w, U = np.linalg.eigh(sum(best[1][i] * P[p][h] for i, p in enumerate(primes))); u = U[:, 0] ** 2; j = np.argsort(-u)[:3]
    print(f"   zoom, {name} half: the falling state is made of near-null states " + ", ".join(f"#{k + 1} of {len(u)} (eigenvalue {EV[h][k]:.1e}, share {u[k]:.2f})" for k in j))
json.dump(dict(c=c, N=N, TOL=TOL, primes=primes, null=[V[0].shape[1], V[1].shape[1]], best=best[0], direction=best[1].tolist(), check=chk, spectrum=spec, directions=[v.tolist() for v in done]), open(f"rigid_{c:g}.json", "w"))
