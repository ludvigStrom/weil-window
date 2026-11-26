"""The bottom of the ordinary spectrum of the odd half, from the zeta zeros alone, over a range of windows (notes/questions.md, Q3).
Q = 2 A A^T as in level_from_zeros.py, modes up to frequency FMAX * 2 pi c, zeros from zeta_zeros_to_6000.json.
Ordinary = eigenvalue above 0.12 (windows are taken at integer + 0.2, where the newest tiny level is near 0.05).
Usage: python3 mu_from_zeros.py FMAX c1 c2 ...   (appends to mu_from_zeros.json)"""
import sys, json, math, os
import numpy as np
FMAX = float(sys.argv[1]); zeros = np.array(json.load(open("zeta_zeros_to_6000.json"))); gaps = np.diff(zeros); out = []
for c in [float(x) for x in sys.argv[2:]]:
    L = math.log(c); W = 2 * math.pi * c; N = int(FMAX * W * L / (2 * math.pi)); w = 2 * math.pi * (np.arange(N) + 1) / L; z = zeros[zeros < 1.6 * w[-1]]
    d = z[None, :] - w[:, None]; s = z[None, :] + w[:, None]
    with np.errstate(divide="ignore", invalid="ignore"): A = math.sqrt(2 / L) * (np.where(np.abs(d) < 1e-12, L / 2, np.sin(d * L / 2) / d) - np.sin(s * L / 2) / s)
    E, V = np.linalg.eigh(2 * A @ A.T); idx = [j for j in range(N) if E[j] > 0.12][:3]; rec = dict(c=c, N=N, states=[])
    for j in idx:
        cw = np.cumsum(V[:, j] ** 2); med = float(w[np.searchsorted(cw, 0.5)]); i = int(np.searchsorted(zeros, med)) - 1
        above = (zeros[:-1] > W) & (zeros[1:] < 0.8 * w[-1]); rank = int((np.where(above, gaps, 0) > gaps[i]).sum()) + 1 if above[i] else 0
        rec["states"].append(dict(eps=float(E[j]), freq=med, z0=float(zeros[i]), z1=float(zeros[i + 1]), rank=rank))
    out.append(rec); s0 = rec["states"][0]
    print(f"c = {c:6.1f}: mu = {s0['eps']:.4f}  (x sqrt(c)/2 = {s0['eps'] * math.sqrt(c) / 2:.2f}) at {s0['freq']:.0f}, gap {s0['z0']:.1f}-{s0['z1']:.1f} = {s0['z1'] - s0['z0']:.2f} = {(s0['z1'] - s0['z0']) * math.log(s0['freq'] / 2 / math.pi) / 2 / math.pi:.2f} mean spacings, rank {s0['rank'] or 'below W'}", flush=True)
R = json.load(open("mu_from_zeros.json")) if os.path.exists("mu_from_zeros.json") else []; json.dump(R + out, open("mu_from_zeros.json", "w"))
