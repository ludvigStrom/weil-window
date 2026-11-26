"""Shared tool: given matrices P_i (derivatives of a form restricted to its near-null space, as a list over halves),
find (a) the exactly flat directions (sum d_i P_i = 0) and (b) the softness spectrum: max over unit d of
lambda_min(sum d_i P_i), then the same orthogonal to what was found. Subgradient ascent on the sphere (the function is concave)."""
import numpy as np
def lam(P, d): return min(np.linalg.eigvalsh(sum(x * Pi[h] for x, Pi in zip(d, P)))[0] for h in range(len(P[0])))
def flat(P, tol=1e-8):
    A = np.array([np.concatenate([Pi[h].ravel() for h in range(len(Pi))]) for Pi in P]).T; U, s, Vt = np.linalg.svd(A, full_matrices=True)
    s = np.concatenate([s, np.zeros(len(P) - len(s))]); return s, Vt[s < tol * s.max()]
def best(P, B, rng, starts=12, iters=3000):
    """max of lam over unit vectors in the column span of B (orthonormal columns)"""
    out = (-1e99, None)
    for t in range(starts):
        y = rng.normal(size=B.shape[1]); y /= np.linalg.norm(y); top = (-1e99, None)
        for k in range(iters):
            d = B @ y; vals = []
            for h in range(len(P[0])):
                w, U = np.linalg.eigh(sum(x * Pi[h] for x, Pi in zip(d, P))); vals.append((w[0], U[:, 0], h))
            w0, u, h = min(vals, key=lambda t: t[0])
            if w0 > top[0]: top = (w0, d.copy())
            gsub = B.T @ np.array([u @ Pi[h] @ u for Pi in P]); gsub -= (gsub @ y) * y; nrm = np.linalg.norm(gsub)
            if nrm < 1e-14: break
            y = y + (0.3 / (1 + k) ** 0.6) * gsub / nrm; y /= np.linalg.norm(y)
        if top[0] > out[0]: out = top
    return out
def spectrum(P, rng, skip=None, **kw):
    n = len(P); done = [] if skip is None else [v for v in skip]; res = []
    while len(done) < n:
        Bc = np.linalg.svd(np.eye(n) - sum(np.outer(v, v) for v in done))[0][:, :n - len(done)] if done else np.eye(n)
        val, d = best(P, Bc, rng, **kw); res.append((val, d)); done.append(d / np.linalg.norm(d))
    return res
