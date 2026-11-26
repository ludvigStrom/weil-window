"""Scan of the image detector of horizon_image.py over depth, from the saved eigenvectors (no form needed).
E(s) = [4 D(s/2) - D(s)] / 3, D(s) = u(w - s) - 2 u(w) + u(w + s), in units of beta, for the states listed.
Usage: python3 horizon_scan.py file.npz s j1,j2,... w_lo w_hi steps"""
import sys, math, numpy as np
Z = np.load(sys.argv[1]); U = Z["U"]; L = float(Z["L"]); beta = Z["beta"]; s = float(sys.argv[2]); js = [int(x) - 1 for x in sys.argv[3].split(",")]
N = U.shape[0]; om = 2 * math.pi * (np.arange(N) + 1) / L
def u(t): return math.sqrt(2 / L) * np.sin(om * (L / 2 - t)) @ U
sgn = np.sign(u(0.004)); c = math.exp(L)
for w in np.linspace(float(sys.argv[4]), float(sys.argv[5]), int(sys.argv[6])):
    D = lambda h: (u(w - h) - 2 * u(w) + u(w + h)) * sgn; E = (4 * D(s / 2) - D(s)) / 3
    print(f"  depth {w:.4f} (c e^-depth = {c * math.exp(-w):.3f}) | " + "  ".join(f"{E[j] / beta[j]:8.3f}" for j in js))
