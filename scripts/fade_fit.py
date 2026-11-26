"""Pinning down the power of p in the fade-in formula (notes/questions.md).
WKB for the prolate equation  -((c - x^2) psi')' + (2 pi)^2 c x^2 psi = chi psi  (c = lambda^2), with chi ~ 2 pi c (2n + 1):
   psi_n(y) ~ [ (4 pi^2 y^2 - 2 pi (2n+1)) (1 - y^2/c) ]^(-1/4) exp( - int_{y_n}^{y} sqrt( (4 pi^2 x^2 - 2 pi (2n+1)) / (1 - x^2/c) ) dx ).
Model: strength_p = (ln p / sqrt p) psi_a(sqrt p) psi_b(sqrt p) x const(window). Which indices (a, b) leave the flattest remainder?
Data: first-zero sensitivities from sens_fast_*.json (primes with one power in the window)."""
import json, math, glob
import numpy as np
from scipy.integrate import quad
rows = []
for fn in sorted(glob.glob("sens_fast_*.json")):
    c = float(fn.split("_")[-1][:-5]); d = json.load(open(fn)); rows += [(c, int(p), s) for p, s in d.items() if int(p) ** 2 > c]
c = np.array([r[0] for r in rows]); p = np.array([r[1] for r in rows], float); S = np.array([r[2] for r in rows]); cs = sorted(set(c)); y0 = S - np.log10(np.log(p) / np.sqrt(p))
def logpsi(n, y, cc):
    k2 = lambda x: (4 * math.pi ** 2 * x * x - 2 * math.pi * (2 * n + 1)) / (1 - x * x / cc); yt = math.sqrt((2 * n + 1) / (2 * math.pi))
    if y <= yt: return None
    return (-quad(lambda x: math.sqrt(max(k2(x), 0)), yt, y)[0] - 0.25 * math.log((4 * math.pi ** 2 * y * y - 2 * math.pi * (2 * n + 1)) * (1 - y * y / cc))) / math.log(10)
def rms(a, b):
    m = [None if (logpsi(a, math.sqrt(pp), cc) is None or logpsi(b, math.sqrt(pp), cc) is None) else logpsi(a, math.sqrt(pp), cc) + logpsi(b, math.sqrt(pp), cc) for pp, cc in zip(p, c)]
    if any(x is None for x in m): return None
    r = y0 - np.array(m); X = np.column_stack([(c == q).astype(float) for q in cs]); th = np.linalg.lstsq(X, r, rcond=None)[0]; return float(np.sqrt(np.mean((X @ th - r) ** 2))), th, r
print(f"{len(rows)} points, windows {cs}; spread of the data {y0.max() - y0.min():.1f} decades")
res = sorted((rms(a, b)[0], a, b) for a in range(0, 13) for b in range(a, 13) if rms(a, b)); print("best index pairs (rms in decades, one free constant per window):")
for r_, a, b in res[:8]: print(f"   n = {a}, {b}: {r_:.3f}")
print("for comparison: n = 0, 0: %.3f;  n = 2, 2: %.3f;  n = 6, 6: %.3f;  n = 2, 6: %.3f" % (rms(0, 0)[0], rms(2, 2)[0], rms(6, 6)[0], rms(2, 6)[0]))
r_, a, b = res[0]; _, th, r = rms(a, b); print(f"constants per window for n = {a}, {b}: " + ", ".join(f"{x:.2f}" for x in th))
for cc in cs: print(f"   c = {cc}: remainders " + ", ".join(f"{int(pp)}: {rr - th[cs.index(cc)]:+.2f}" for pp, q, rr in zip(p, c, r) if q == cc))
