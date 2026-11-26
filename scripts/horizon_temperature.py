"""Test 2 of notes/picture.md round 5: is the ladder of tiny eigenvalues a Boltzmann law at the horizon's temperature?
Fermi gas in a box, W = 2 pi c, prolate eigenvalue chi_n, q = chi_n / W^2.
  Horizon law (Bisognano-Wichmann / Hawking, temperature 1/(2 pi) in boost units; the entanglement Hamiltonian is the prolate
  operator rescaled):                 ln(1/leak) = (pi / 2W) (W^2 - chi_n)             [linear in the prolate energy]
  Full WKB (notes/margin.md part 12):       ln(1/leak) = 2 W T(q),  T = E(k') - q K(k')      [reduces to the line above as q -> 1]
  Fuchs, fixed n, large c:            delta_n = 2 sqrt(pi) 8^n (2 pi)^(n+1/2) / n! * c^(n+1/2) e^(-4 pi c)
Odd Weil states use n = 6, 10, 14, ... (Connes-Consani). Printed: ln(1/eps_j) measured against the three, and the step between
neighbours against 4 pi K(k')/K(k) (WKB) and against the thermal step (pi/2W) (chi_{n+4} - chi_n)."""
import json, math
from scipy.special import pro_cv, ellipk, ellipe
eps = {7.021032: [1.2599e-24, 1.0567e-18, 2.6870e-13, 6.9523e-09, 9.5656e-05, 7.3626e-02]}
for c in (11.011, 13.013): eps[c] = [e for e in json.load(open(f"edge_rates_{c}.json"))["eps"] if e < 0.2]
for c, E in eps.items():
    W = 2 * math.pi * c; print(f"c = {c}, W = {W:.2f}"); print("   j | n | q | ln(1/eps) measured | horizon law | WKB | Fuchs || step to next: measured | thermal step | WKB step 4 pi K'/K")
    rows = []
    for j, e in enumerate(E):
        n = 6 + 4 * j; chi = pro_cv(0, n, W); q = chi / W ** 2
        hor = math.pi / (2 * W) * (W * W - chi); wkb = 2 * W * (ellipe(1 - q) - q * ellipk(1 - q)) if 0 < q < 1 else float("nan")
        fu = -(math.log(2 * math.sqrt(math.pi)) + n * math.log(8) + (n + .5) * math.log(2 * math.pi) - math.lgamma(n + 1) + (n + .5) * math.log(c) - 4 * math.pi * c)
        rows.append((n, q, -math.log(e), hor, wkb, fu, chi))
    for j, (n, q, m, hor, wkb, fu, chi) in enumerate(rows):
        if j + 1 < len(rows):
            n2, q2, m2, hor2, wkb2, fu2, chi2 = rows[j + 1]; qm = (q + q2) / 2
            step = f"{m - m2:6.2f} | {hor - hor2:6.2f} | {4 * math.pi * ellipk(1 - qm) / ellipk(qm) if qm < 1 else float('nan'):6.2f}"
        else: step = ""
        print(f"   {j + 1:2d} | {n:2d} | {q:.3f} | {m:7.2f} | {hor:7.2f} | {wkb:7.2f} | {fu:7.2f} || {step}")
