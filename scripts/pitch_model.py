"""Chain model, second attempt: sites labelled by frequency (notes/questions.md).
Sites = eigenvectors of the frequency operator (diag of the mode frequencies) compressed to the near-null space V,
each half separately. On the whole line a shift by a is multiplication by 2 cos(omega a) in frequency. Test stated
beforehand: >= 80% of each prime power's action on the diagonal, and diagonal values within 10% of 2 cos(nu_j a).
Also: where do the sites sit relative to the zeta zeros?   Usage: python3 pitch_model.py N c   (needs rigid_P_<c>_<N>.npz)"""
import sys, math, os, json
import numpy as np
N = int(sys.argv[1]); c = float(sys.argv[2]); TOL = float(os.environ.get("TOL", 1e-2)); L = math.log(c); W = 2 * math.pi * c
Z = dict(np.load(f"rigid_P_{c:g}_{N}.npz")); m = np.arange(N + 1); om = 2 * math.pi * m / L; xg, wg = np.polynomial.legendre.leggauss(3000)
def basis(u):
    C = np.cos(np.outer(u, om)) * np.where(m == 0, 1 / math.sqrt(L), math.sqrt(2 / L)); S = np.sin(np.outer(u, om[1:])) * math.sqrt(2 / L); return C, S
def A(a):
    lo, hi = -L / 2, L / 2 - a; u = (hi - lo) / 2 * xg + (hi + lo) / 2; wq = (hi - lo) / 2 * wg; C1, S1 = basis(u); C2, S2 = basis(u + a)
    Ge = (C1 * wq[:, None]).T @ C2; Go = (S1 * wq[:, None]).T @ S2; return Ge + Ge.T, Go + Go.T
zeros = np.array(json.load(open("zeta_zeros_to_2300.json"))); pps = sorted((q, p) for p in range(2, int(c) + 1) if all(p % d for d in range(2, int(p ** 0.5) + 1)) for q in [p ** k for k in range(1, 12)] if q < c)
print(f"c = {c}: 2 pi c = {W:.1f}; zeros below it: {int((zeros < W).sum())}")
for h, name, Q0, freq in ((0, "even", Z["Q0e"], om), (1, "odd", Z["Q0o"], om[1:])):
    E, U = np.linalg.eigh(Q0); V = U[:, E < TOL]; nu, R = np.linalg.eigh(V.T @ np.diag(freq) @ V); S = V @ R
    width = [math.sqrt(max(s @ (freq ** 2 * s) - n * n, 0)) for s, n in zip(S.T, nu)]
    print(f"   {name} half, {len(nu)} sites at frequencies: " + ", ".join(f"{x:.1f}" for x in nu)); print(f"      width of a site in frequency: {min(width):.1f} to {max(width):.1f}; spacing {np.diff(nu).min():.1f} to {np.diff(nu).max():.1f}")
    nearest = [float(np.min(np.abs(zeros - x))) for x in nu]; sp = [float(np.diff(zeros)[max(np.searchsorted(zeros, x) - 1, 0)]) for x in nu]
    print("      distance from each site to the nearest zero, as a fraction of the local gap between zeros (0.5 = mid-gap): " + ", ".join(f"{d / s:.2f}" for d, s in zip(nearest, sp)))
    arch = V.T @ Q0 @ V
    for q, p in pps:
        a = math.log(q); M = S.T @ A(a)[h] @ S; d = np.diag(M); share = (d ** 2).sum() / (M ** 2).sum(); pred = 2 * np.cos(nu * a) * max(0.0, 1 - a / L); err = np.abs(d - pred).mean()
        arch = arch + math.log(p) / math.sqrt(q) * (V.T @ A(a)[h] @ V); print(f"      {q:3d}: share on-site {share:.2f}; on-site values against 2 (1 - a/L) cos(nu a): mean abs difference {err:.2f} (values range {d.min():+.2f} to {d.max():+.2f})")
    Ms = R.T @ arch @ R; print(f"      background: {(np.diag(Ms) ** 2).sum() / (Ms ** 2).sum():.2f} on-site; values " + ", ".join(f"{x:.2f}" for x in np.diag(Ms)))
