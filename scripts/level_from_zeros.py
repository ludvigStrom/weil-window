"""Question 3 (notes/questions.md): predict the low ordinary levels of the odd half from the zeta zeros alone.
Weil's form is the sum over zeros: for odd real f on the window [-L/2, L/2], Q(f) = 2 sum_{gamma > 0} |int f(x) sin(gamma x) dx|^2.
In the basis s_k = sqrt(2/L) sin(w_k x), w_k = 2 pi k / L:   Q = 2 A A^T,  A[k, gamma] = int s_k(x) sin(gamma x) dx  (closed form).
(a) Whole basis, all zeros below 2300: compare the ordinary eigenvalues with those of the true form (flat_track.json).
(b) Local model: only the modes within BAND of a gap's centre and the zeros within REACH of it: lowest eigenvalue.
Usage: python3 level_from_zeros.py c N [BAND REACH]"""
import sys, json, math
import numpy as np
c = float(sys.argv[1]); N = int(sys.argv[2]); BAND = float(sys.argv[3]) if len(sys.argv) > 3 else 25.0; REACH = float(sys.argv[4]) if len(sys.argv) > 4 else 120.0
L = math.log(c); W = 2 * math.pi * c; zeros = np.array(json.load(open("zeta_zeros_to_2300.json"))); w = 2 * math.pi * (np.arange(N) + 1) / L
def A(wk, gam):
    d = gam[None, :] - wk[:, None]; s = gam[None, :] + wk[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        t1 = np.where(np.abs(d) < 1e-12, L / 2, np.sin(d * L / 2) / d)
    return math.sqrt(2 / L) * (t1 - np.sin(s * L / 2) / s)
true = min((r for r in json.load(open("flat_track.json")) if abs(r["c"] - c) < 1e-9), key=lambda r: r["ntiny"])
tv = [(s["eps"], s["band"][1]) for s in true["states"] if not (s["eps"] < 0.1 and s["edge"] > 0.02)]
M = A(w, zeros); Q = 2 * M @ M.T; E, V = np.linalg.eigh(Q); cw = np.cumsum(V ** 2, axis=0)
print(f"c = {c}, W = {W:.0f}, basis to frequency {w[-1]:.0f}, {len(zeros)} zeros below 2300")
print("(a) ordinary levels: true form (ball arithmetic) | from the zeros alone, with the state's median frequency")
ev = [(E[j], w[np.searchsorted(cw[:, j], 0.5)]) for j in range(N) if E[j] > 0.1][:8]
for (a, fa), (b, fb) in zip(tv[:8], ev): print(f"     {a:.4f} at {fa:6.1f}   |   {b:.4f} at {fb:6.1f}")
print(f"(b) local model: modes within {BAND:g} of the gap centre, zeros within {REACH:g}")
gaps = np.diff(zeros); usable = (zeros[:-1] > 0.9 * W) & (zeros[1:] < 0.85 * w[-1]); order = np.argsort(-np.where(usable, gaps, 0))[:8]
print("     gap (zeros) | width | local-model level | true states sitting there")
for i in order:
    mid = (zeros[i] + zeros[i + 1]) / 2; km = np.abs(w - mid) < BAND; zm = np.abs(zeros - mid) < REACH
    Ml = A(w[km], zeros[zm]); lam = np.linalg.eigvalsh(2 * Ml @ Ml.T)[0]
    here = [f"{a:.4f}" for a, fa in tv if zeros[i] - 0.5 <= fa <= zeros[i + 1] + 0.5]
    print(f"     {zeros[i]:8.2f} - {zeros[i + 1]:8.2f} | {gaps[i]:.2f} | {lam:.4f} | {', '.join(here) if here else '-'}")
