"""What are the bumps in the level's rate of fall? (notes/questions.md)
Data: level_curve.json (fixed N). Tests: (1) fold the rate by the fractional part of c (new near-null states arrive once
per unit of c: odd half at integers, even half at half-integers); (2) periodogram of the rate in several variables
(c, sqrt c, ln c, number of zeros below 2 pi c) to see in which one, if any, the bumps are periodic.
Usage: python3 bumps.py [N]"""
import json, math, sys
import numpy as np
NN = sys.argv[1] if len(sys.argv) > 1 else "80"
d = json.load(open("level_curve.json")); pts = sorted((float(k.split("|")[0]), v) for k, v in d.items() if k.endswith("|" + NN)); c = np.array([p[0] for p in pts]); e = np.array([p[1][0] for p in pts]); o = np.array([p[1][1] for p in pts])
se = -np.diff(e) / np.diff(c); so = -np.diff(o) / np.diff(c); cm = (c[1:] + c[:-1]) / 2; print(f"{len(cm)} intervals, c from {c[0]} to {c[-1]}")
def detrend(y): p = np.polyfit(cm, y, 1); return y - np.polyval(p, cm)
de, do = detrend(se), detrend(so)
print("(1) folded by the fractional part of c: mean of the detrended rate in five bins (standard error about %.2f)" % (de.std() / math.sqrt(len(cm) / 5)))
fr = cm % 1
for name, y in (("even", de), ("odd", do)):
    print(f"   {name}: " + "  ".join(f"[{a:.1f},{a + .2:.1f}): {y[(fr >= a) & (fr < a + .2)].mean():+.3f}" for a in (0, .2, .4, .6, .8)) + f"    variance explained by the fold: {1 - np.mean([(y[(fr >= a) & (fr < a + .2)] - y[(fr >= a) & (fr < a + .2)].mean()).var() for a in (0, .2, .4, .6, .8)]) / y.var():.2f}")
print(f"   correlation between the even and odd rates: {np.corrcoef(de, do)[0, 1]:+.2f}; with the odd one shifted by half a unit: {np.corrcoef(de[10:], do[:-10])[0, 1]:+.2f}; by a quarter: {np.corrcoef(de[5:], do[:-5])[0, 1]:+.2f}")
print("(2) periodogram: the three strongest periods, and the share of the variance in the strongest, in each variable")
def pgram(x, y, name):
    span = x[-1] - x[0]; fs = np.linspace(1 / span, len(x) / (2.2 * span), 600); P = np.array([abs(np.sum((y - y.mean()) * np.exp(-2j * np.pi * f * x))) ** 2 for f in fs]); P /= (len(x) * (y - y.mean()).var() * len(x) / 4)
    idx = [i for i in range(1, len(fs) - 1) if P[i] > P[i - 1] and P[i] > P[i + 1]]; idx = sorted(idx, key=lambda i: -P[i])[:3]
    return f"   {name:34s} " + "; ".join(f"period {1 / fs[i]:.3f} (share {P[i]:.2f})" for i in idx)
Nz = cm * np.log(cm) - cm
for nm, y in (("even", de), ("odd", do)):
    print(f"  {nm} half:")
    for name, x in (("c", cm), ("sqrt c = lambda", np.sqrt(cm)), ("ln c = L", np.log(cm)), ("zeros below 2 pi c: c ln c - c", Nz), ("c^2", cm ** 2)): print(pgram(x, y, name))
