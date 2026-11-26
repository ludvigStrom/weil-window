"""Is the form, on its near-null states, a hopping model? (notes/questions.md, chain model)
Sites = eigenvectors of the position operator u compressed to the near-null space V (even and odd halves together):
the most localised basis V allows. In that basis: is each prime's operator A_a concentrated on pairs of sites a apart?
Is the background (form + all prime terms) concentrated on the diagonal?
Success test stated beforehand: >= 80% of each prime's action on links of the right length; background mostly on-site.
Usage: python3 chain_model.py N c      (needs rigid_P_<c>_<N>.npz; env TOL default 1e-2)"""
import sys, math, os
import numpy as np
N = int(sys.argv[1]); c = float(sys.argv[2]); TOL = float(os.environ.get("TOL", 1e-2)); L = math.log(c)
Z = dict(np.load(f"rigid_P_{c:g}_{N}.npz")); m = np.arange(N + 1); om = 2 * math.pi * m / L; xg, wg = np.polynomial.legendre.leggauss(3000)
def basis(u):
    C = np.cos(np.outer(u, om)) * np.where(m == 0, 1 / math.sqrt(L), math.sqrt(2 / L)); S = np.sin(np.outer(u, om[1:])) * math.sqrt(2 / L); return np.hstack([C, S])
def A(a):
    lo, hi = -L / 2, L / 2 - a; u = (hi - lo) / 2 * xg + (hi + lo) / 2; wq = (hi - lo) / 2 * wg; G = (basis(u) * wq[:, None]).T @ basis(u + a); return G + G.T
u = L / 2 * xg; B = basis(u); X = (B * (L / 2 * wg * u)[:, None]).T @ B                         # position operator in the full basis
Ve = np.linalg.eigh(Z["Q0e"]); Vo = np.linalg.eigh(Z["Q0o"]); ne, no = int((Ve[0] < TOL).sum()), int((Vo[0] < TOL).sum())
V = np.zeros((2 * N + 1, ne + no)); V[:N + 1, :ne] = Ve[1][:, :ne]; V[N + 1:, ne:] = Vo[1][:, :no]
Q0 = np.zeros((2 * N + 1, 2 * N + 1)); Q0[:N + 1, :N + 1] = Z["Q0e"]; Q0[N + 1:, N + 1:] = Z["Q0o"]
pos, R = np.linalg.eigh(V.T @ X @ V); S = V @ R                                                  # sites: columns of S, at positions pos
width = [math.sqrt(max(s @ (B * (L / 2 * wg * u * u)[:, None]).T @ B @ s - p * p, 0)) for s, p in zip(S.T, pos)]
gap = np.diff(pos); print(f"c = {c}: window [-{L / 2:.3f}, {L / 2:.3f}], {len(pos)} sites")
print("   site positions: " + ", ".join(f"{p:+.2f}" for p in pos)); print(f"   spacing between sites: {gap.min():.3f} to {gap.max():.3f} (median {np.median(gap):.3f}); width of a site's wave: {min(width):.3f} to {max(width):.3f}")
pps = sorted((q, p) for p in range(2, int(c) + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1)) for q in [p ** k for k in range(1, 12)] if q < c)
dist = np.abs(pos[:, None] - pos[None, :]); tol = float(np.median(gap)); rng = np.random.default_rng(0); O = np.linalg.qr(rng.normal(size=(len(pos), len(pos))))[0]; arch = V.T @ Q0 @ V
print(f"   share of each prime power's action on pairs of sites within {tol:.3f} of the right distance (and the same in a random basis of V):")
for q, p in pps:
    a = math.log(q); M = S.T @ A(a) @ S; mask = np.abs(dist - a) < tol; share = (M[mask] ** 2).sum() / (M ** 2).sum(); Mr = O.T @ M @ O; arch = arch + math.log(p) / math.sqrt(q) * (V.T @ A(a) @ V)
    print(f"      {q:3d} (a = {a:.3f}): {share:.2f}   (random basis: {(Mr[mask] ** 2).sum() / (Mr ** 2).sum():.2f}; share of site pairs that qualify: {mask.mean():.2f})")
Ms = R.T @ arch @ R; print(f"   background (form + all prime terms) in the site basis: {(np.diag(Ms) ** 2).sum() / (Ms ** 2).sum():.2f} of its weight on the diagonal, {sum((np.diag(Ms, k) ** 2).sum() for k in (-1, 1)) / (Ms ** 2).sum():.2f} on nearest neighbours; on-site values " + ", ".join(f"{x:.2f}" for x in np.diag(Ms)))
