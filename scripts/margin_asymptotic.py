"""The margin criterion for large n, on paper (notes/margin.md, part 12).
Prolate WKB, W = 2 pi c, q = chi_n / W^2:   index  n + 1/2 = (2W/pi) X(q),  X = E(k) - (1 - q) K(k),  k^2 = q;
leak 1 - lambda_n ~ exp(-2 W T(q)),  T = E(k') - q K(k'),  k'^2 = 1 - q.   Rate at fixed n, with dW/da = 2W and Legendre's relation:
    -d ln(leak)/da = 2 pi W / K(k) = 4 pi^2 c / K(sqrt q),      so beta^2/eps = 2 pi^2 c / K.
Sum over the odd states (every fourth mode, c of them): S = 2 pi^2 c^2 int_0^1 dx / K = pi^2 c^2   (dx = (K/2) dq).
With C = W sqrt(2 (1 - q)), w = mu / c^2 and the Bessel factor 2 I1(z)/z, the sum is (I0(a) - 1)/4 with a = 4 pi sqrt(mu):
    2 lambda [ floor + (I0(4 pi sqrt(mu)) - 1) / (4 l) ] = 1,   l = w / I_kappa(w)  (about ln(1/(kappa w)))."""
import json, math, numpy as np
from scipy.special import ellipk, i0, pro_cv
from scipy.optimize import brentq
KAP = 2 * math.pi * math.exp(0.5772156649); FLOOR = 0.2
had = {7.021032: [79.51, 79.98, 72.73, 70.97, 60.89, 51.31]}
print("1. Rates beta^2/eps: measured (deepest first) against 2 pi^2 c / K(sqrt q_n), n = 10, 14, 18, ...")
for c in (7.021032, 11.011, 13.013):
    if c in had: S = had[c]
    else: R = json.load(open(f"edge_rates_{c}.json")); S = [b / e for b, e in zip(R["beta2"], R["eps"]) if e < 0.2]
    W = 2 * math.pi * c; pred = []
    for j in range(len(S)):
        q = pro_cv(0, 10 + 4 * j, W) / W ** 2; pred.append(2 * math.pi ** 2 * c / ellipk(q) if q < 1 else float("nan"))
    print(f"   c = {c}: measured " + " ".join(f"{s:.0f}" for s in S)); print(f"   {'':12s} formula  " + " ".join(f"{p:.0f}" for p in pred) + f"     sum measured {sum(S):.0f}, pi^2 c^2 = {math.pi ** 2 * c * c:.0f}")
def Ik(w):
    t = (np.arange(40000) + 0.5) / 40000 * w; return float((1 / np.sqrt(np.log(1 / (KAP * t)) * np.log(1 / (KAP * (w - t))))).sum() * w / 40000)
print("2. Margin from the closed criterion:  n | mu = w* n^2 | w* | (model of part 11, measured)")
ref = {3: (0.0211, 0.0163), 5: (0.0055, 0.0051), 7: (0.0026, 0.0025), 11: (0.0010, None), 13: (0.0007, None)}
for n in (3, 5, 7, 11, 13, 101, 1009, 10007, 100003, 1000003):
    lam = math.log(n) / math.sqrt(n)
    def g(mu): w = mu / n ** 2; return FLOOR + (i0(4 * math.pi * math.sqrt(mu)) - 1) * Ik(w) / (4 * w) - 1 / (2 * lam)
    mu = brentq(g, 1e-4, min(5.0, 0.5 * n * n / KAP)); print(f"   {n:8d} | {mu:.4f} | {mu / n ** 2:.3e} | {ref.get(n, '')}")
