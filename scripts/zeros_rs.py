"""Zeta zeros up to height T from the Riemann-Siegel formula in double precision (one correction term).
Good to about 0.005 in the ordinate, which is enough to tell which gap a state sits in. Writes zeta_zeros_to_<T>.json.
Usage: python3 zeros_rs.py [T]     (default 2300)"""
import sys, json, math
import numpy as np
T = float(sys.argv[1]) if len(sys.argv) > 1 else 2300.0
def Z(t):
    a = np.sqrt(t / (2 * np.pi)); n_max = np.floor(a).astype(int); p = a - n_max
    theta = t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)
    out = np.zeros_like(t)
    for n in range(1, int(n_max.max()) + 1):
        out += np.where(n <= n_max, 2 * np.cos(theta - t * np.log(n)) / np.sqrt(n), 0)
    return out + (-1.0) ** (n_max - 1) * (t / (2 * np.pi)) ** -0.25 * np.cos(2 * np.pi * (p * p - p - 1 / 16)) / np.cos(2 * np.pi * p)
ts = np.arange(10, T, 0.01); v = Z(ts); i = np.where(v[:-1] * v[1:] < 0)[0]
zeros = (ts[i] - v[i] * (ts[i + 1] - ts[i]) / (v[i + 1] - v[i])).tolist()
expected = T / (2 * math.pi) * math.log(T / (2 * math.pi * math.e)) + 7 / 8
print(f"{len(zeros)} zeros below {T:g}; the counting formula gives about {expected:.1f}")
json.dump(zeros, open(f"zeta_zeros_to_{T:g}.json", "w"))
