"""Why is a prime held as hard as it is? (notes/questions.md, pattern stiffness)
The term of a prime power at log-position a is  -w * A_a,  A_a(f, g) = int f(u) g(u+a) + g(u) f(u+a) du  over the window
[-L/2, L/2] (only the stretch where both points are inside). So
   weight tent   = w * extreme eigenvalues of A_a        on the near-null space V,
   position tent = w * extreme eigenvalues of dA_a/da    on V   (for a prime with one power in the window).
This script (1) checks that identification against the matrices from rigid.py / rigid_weights.py, and
(2) prints the two profiles as functions of a: how hard a mirror placed anywhere in the window would be held.
Usage: python3 stiff_profile.py N c      (needs rigid_P_<c>_<N>.npz and rigid_W_<c>_<N>.npz; env TOL default 1e-2)"""
import sys, math, os, json
import numpy as np
N = int(sys.argv[1]); c = float(sys.argv[2]); TOL = float(os.environ.get("TOL", 1e-2)); L = math.log(c); W = 2 * math.pi * c
Z = dict(np.load(f"rigid_P_{c:g}_{N}.npz")); Wt = dict(np.load(f"rigid_W_{c:g}_{N}.npz")) if os.path.exists(f"rigid_W_{c:g}_{N}.npz") else {}; m = np.arange(N + 1); om = 2 * math.pi * m / L
xg, wg = np.polynomial.legendre.leggauss(3000)
def basis(u):                                    # rows: points; even block (N+1 columns), odd block (N columns)
    C = np.cos(np.outer(u, om)) * np.where(m == 0, 1 / math.sqrt(L), math.sqrt(2 / L)); S = np.sin(np.outer(u, om[1:])) * math.sqrt(2 / L); return C, S
def A(a):
    lo, hi = -L / 2, L / 2 - a; u = (hi - lo) / 2 * xg + (hi + lo) / 2; wq = (hi - lo) / 2 * wg; C1, S1 = basis(u); C2, S2 = basis(u + a)
    Ge = (C1 * wq[:, None]).T @ C2; Go = (S1 * wq[:, None]).T @ S2; return Ge + Ge.T, Go + Go.T
dA = lambda a, h=1e-5: [(x - y) / (2 * h) for x, y in zip(A(a + h), A(a - h))]
V = []
for Mx in (Z["Q0e"], Z["Q0o"]):
    E, U = np.linalg.eigh(Mx); V.append(U[:, E < TOL])
ext = lambda Ms: (min(np.linalg.eigvalsh(V[h].T @ Ms[h] @ V[h])[0] for h in range(2)), max(np.linalg.eigvalsh(V[h].T @ Ms[h] @ V[h])[-1] for h in range(2)))
print(f"c = {c}: L = {L:.3f}, 2 pi c = {W:.1f}; near-null space {V[0].shape[1]} + {V[1].shape[1]}")
# (1) the identification, on a prime with a single power in the window
for p in [q for q in (5, 7, 11, 13, 17, 19, 23) if q <= c and q * q > c and f"De{q}" in Z]:
    a = math.log(p); w = a / math.sqrt(p); Aa = A(a); r = [np.linalg.norm(Wt[k] + w * Aa[h]) / np.linalg.norm(Wt[k]) for h, k in ((0, f"We{p}"), (1, f"Wo{p}"))] if Wt else [float("nan")]
    dAa = dA(a); model = [-w * (dAa[h] + (1 / a - 0.5) * Aa[h]) for h in range(2)]          # d/d(delta): position a = ln p + delta, weight ln x / sqrt x
    r2 = [np.linalg.norm(V[h].T @ (Z[k] - model[h]) @ V[h]) / np.linalg.norm(V[h].T @ Z[k] @ V[h]) for h, k in ((0, f"De{p}"), (1, f"Do{p}"))]
    eA, eD = ext(Aa), ext(dAa); print(f"   prime {p}: weight matrix = -w A_a to {max(r):.1e}; position matrix = -w (dA/da + (1/a - 1/2) A) to {max(r2):.1e} on V;  extremes of A_a on V: {eA[0]:+.3f}, {eA[1]:+.3f};  of dA/da: {eD[0]:+.2f}, {eD[1]:+.2f};  w x these = {w * eD[0]:+.1f}, {w * eD[1]:+.1f}")
# (2) the profiles
print("   a/L    a=ln x: x    extremes of A_a on V     extremes of dA_a/da on V    (max |dA/da|) / (2 pi c)")
prof = []
for t in np.linspace(0.05, 0.95, 19):
    a = t * L; eA, eD = ext(A(a)), ext(dA(a)); prof.append(dict(t=float(t), x=math.exp(a), A=eA, dA=eD)); print(f"   {t:.2f}   {math.exp(a):7.3f}      {eA[0]:+.3f}, {eA[1]:+.3f}          {eD[0]:+8.2f}, {eD[1]:+8.2f}            {max(abs(eD[0]), abs(eD[1])) / W:.3f}", flush=True)
json.dump(prof, open(f"stiff_profile_{c:g}.json", "w"))
