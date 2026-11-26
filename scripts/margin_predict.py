"""Margin of a withheld prime from the model of part 10 (notes/margin.md, part 11).
2 lambda [ floor + sum_j (beta_j^2/eps_j) int_0^w P_j(t) P_j(w - t) dt ] = 1,  P_j = I0(C_j sqrt t)/sqrt(ln(1/(kappa t))),
over the near-null odd states (eps < 0.2). beta^2/eps: Hadamard formula at N = 240; C_j: fitted to the profile at N = 360,
or the prolate value sqrt(2 (W^2 - chi_n)), W = 2 pi c, n = 10, 14, 18, ... from the deepest state."""
import json, math, os, numpy as np
from scipy.special import i0, pro_cv
from scipy.optimize import minimize_scalar, brentq
KAP = 2 * math.pi * math.exp(0.5772156649); FLOOR = 0.2
had = {3.049: [2.7877e-04, 1.1813], 5.0255: [4.3156e-13, 7.2900e-08, 1.8806e-03, 2.4264], 7.021032: [1.0017e-22, 8.4516e-17, 1.9544e-11, 4.9339e-07, 5.8249e-03, 3.7780]}
meas = {3: 0.0163, 5: 0.0051, 7: 0.0025}
def rho(w, S, Cs):
    tt = (np.arange(40000) + 0.5) / 40000 * w; lw = 1 / np.sqrt(np.log(1 / (KAP * tt)) * np.log(1 / (KAP * (w - tt))))
    return FLOOR + sum(s * (lw * i0(C * np.sqrt(tt)) * i0(C * np.sqrt(w - tt))).sum() * w / 40000 for s, C in zip(S, Cs))
print("  n | near-null states | S = sum beta^2/eps | C fitted (deepest first) | C prolate | w* (fitted C) | w* (prolate C) | measured | w* n^2 | mean Bessel factor")
for n, c in ((3, 3.049), (5, 5.0255), (7, 7.021032), (11, 11.011), (13, 13.013)):
    D = json.load(open(f"deep_profile_{c}_360.json")); t = np.array(D["t"]); u = np.abs(np.array(D["u"])); eps = D["eps"]
    b2 = had[c] if c in had else json.load(open(f"edge_rates_{c}.json"))["beta2"]; K = sum(e < 0.2 for e in eps); S = [b2[j] / eps[j] for j in range(K)]; Cs = []
    for j in range(K):
        m = (t > 0.003) & (t < (0.06 if j < K - 1 else 0.02)); y = np.log(u[m, j] * np.sqrt(np.log(1 / (KAP * t[m]))))
        def cost(C): r = y - np.log(i0(C * np.sqrt(t[m]))); return ((r - r.mean()) ** 2).sum()
        Cs.append(minimize_scalar(cost, bounds=(0.01, 200), method="bounded").x)
    W = 2 * math.pi * c; Cp = [math.sqrt(max(2 * (W * W - pro_cv(0, 10 + 4 * j, W)), 0.0)) for j in range(K)]
    lam = math.log(n) / math.sqrt(n); need = 1 / (2 * lam); wf = brentq(lambda w: rho(w, S, Cs) - need, 1e-5, 0.05); wp = brentq(lambda w: rho(w, S, Cp) - need, 1e-5, 0.05)
    bes = (rho(wf, S, Cs) - FLOOR) / (rho(wf, S, [0.0] * K) - FLOOR)
    print(f"  {n} | {K} | {sum(S):.0f} | {' '.join(f'{C:.0f}' for C in Cs)} | {' '.join(f'{C:.0f}' for C in Cp)} | {wf:.5f} | {wp:.5f} | {meas.get(n, float('nan'))} | {wf * n * n:.3f} | {bes:.2f}")
