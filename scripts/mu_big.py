"""Bottom of the ordinary spectrum of the odd half from the zeta zeros alone, for large windows (notes/questions.md, Q3 continued).
Same Q = 2 A A^T as mu_from_zeros.py, but only the sine modes with frequency between FLO and FHI times 2 pi c are kept
(FLO = 0 is the full basis; a band is a Rayleigh-Ritz upper bound, checked against the full basis where both fit in memory),
and Q is accumulated over chunks of zeros (all zeros below ZHI * 2 pi c, from zeta_zeros_to_21000.json).
Prints every eigenvalue between 0.02 and EMAX with the state's median frequency and the gap between zeros it sits in.
Usage: python3 mu_big.py FLO FHI ZHI c1 c2 ...   (appends to mu_big.json; env EMAX default 0.3)"""
import sys, json, math, os, time
import numpy as np
FLO, FHI, ZHI = (float(x) for x in sys.argv[1:4]); EMAX = float(os.environ.get("EMAX", 0.3))
zeros = np.array(json.load(open("zeta_zeros_to_21000.json"))); out = []
for c in [float(x) for x in sys.argv[4:]]:
    t0 = time.time(); L = math.log(c); W = 2 * math.pi * c
    k = np.arange(max(1, int(FLO * W * L / (2 * math.pi))), int(FHI * W * L / (2 * math.pi)) + 1); w = 2 * math.pi * k / L; N = len(w)
    z = zeros[zeros < ZHI * W]; assert ZHI * W < zeros[-1], "need more zeros"; Q = np.zeros((N, N))
    for j in range(0, len(z), 2000):
        zz = z[j:j + 2000]; d = zz[None, :] - w[:, None]; s = zz[None, :] + w[:, None]
        with np.errstate(divide="ignore", invalid="ignore"): A = math.sqrt(2 / L) * (np.where(np.abs(d) < 1e-12, L / 2, np.sin(d * L / 2) / d) - np.sin(s * L / 2) / s)
        Q += 2 * A @ A.T
    E, V = np.linalg.eigh(Q); rec = dict(c=c, FLO=FLO, FHI=FHI, ZHI=ZHI, N=N, states=[])
    print(f"c = {c}: W = {W:.0f}, {N} modes, {len(z)} zeros, {int((E < 0.02).sum())} eigenvalues below 0.02, {time.time() - t0:.0f} s", flush=True)
    for j in np.where((E > 0.02) & (E < EMAX))[0][:14]:
        cw = np.cumsum(V[:, j] ** 2); med = float(w[np.searchsorted(cw, 0.5)]); lo, hi = float(w[np.searchsorted(cw, 0.1)]), float(w[np.searchsorted(cw, 0.9)]); i = int(np.searchsorted(zeros, med)) - 1
        g = (zeros[i + 1] - zeros[i]) * math.log(med / 2 / math.pi) / 2 / math.pi
        rec["states"].append(dict(eps=float(E[j]), freq=med, lo=lo, hi=hi, z0=float(zeros[i]), z1=float(zeros[i + 1]), gap=g))
        print(f"   {E[j]:.4f}  (x sqrt(c)/2 = {E[j] * math.sqrt(c) / 2:.2f})  freq {med:7.1f} [{lo:.0f}, {hi:.0f}] = {med / W:.2f} W   gap {zeros[i]:.2f}-{zeros[i + 1]:.2f} = {g:.2f} mean spacings", flush=True)
    out.append(rec)
R = json.load(open("mu_big.json")) if os.path.exists("mu_big.json") else []; json.dump(R + out, open("mu_big.json", "w"))
