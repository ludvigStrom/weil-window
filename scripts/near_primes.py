"""Near primes from Riemann's own function (notes/questions.md, 'Near primes').
As the window grows the vacuum tends to Phi(u) = 2 sum_n (2 pi^2 n^4 e^(9u/2) - 3 pi n^2 e^(5u/2)) exp(-pi n^2 e^(2u)),
the even function whose cosine transform is Xi(s) = xi(1/2 + i s). Then with no window and no fitted number:
   R(a)   = int Phi(u) Phi(u + a) du / int Phi^2            (self-overlap; its transform is Xi^2)
   P_n(a) = (1 / g_n) int_a^inf R(y) sin(g_n (y - a)) dy     (the same with zero g_n divided out)
   rho_n(a) = 1 - (P_n(a) / P_n(0)) / R(a)                   (response of zero n to a mirror at lag a, in units of K_n x level shift)
Compared with the measured ratios at c = 28.2 (response_curve.py KCHK).   Double precision on a grid; no window, no fitted number.   Usage: python3 near_primes.py"""
import numpy as np
du = 5e-4; u = np.arange(-4, 4 + du / 2, du); n = np.arange(1, 7)[:, None]
def Phi(x):
    x = np.abs(x); return 2 * np.sum((2 * np.pi ** 2 * n ** 4 * np.exp(4.5 * x) - 3 * np.pi * n ** 2 * np.exp(2.5 * x)) * np.exp(-np.pi * n ** 2 * np.exp(2 * x)), axis=0)
F = Phi(u); N0 = np.sum(F * F) * du; M = int(round(2.4 / du)); Rg = np.array([np.sum(F[:len(F) - k] * F[k:]) * du / N0 for k in range(M + 1)]); yg = np.arange(M + 1) * du      # R on a grid of lags
def P(g, k0): return np.sum(Rg[k0:] * np.sin(g * (yg[k0:] - yg[k0]))) * du / g
zs = [14.134725141734695, 21.022039638771556, 30.424876125859512]; meas = {1.5: (0.111, 0.044, 0.020), 2: (0.298, 0.139, 0.066), 3: (0.589, 0.367, 0.204), 5: (0.837, 0.686, 0.498), 7: (0.919, 0.831, 0.695)}
P0 = [P(g, 0) for g in zs]; print("check: transform of R at the zeros (should vanish): " + ", ".join(f"{2 * np.sum(Rg * np.cos(g * yg)) * du - du:.1e}" for g in zs) + ";  g^2 P(0) = " + ", ".join(f"{g * g * p:.4f}" for g, p in zip(zs, P0)))
print("   q     R(ln q) (log10)    rho for zeros 1, 2, 4 from Phi alone        measured at c = 28.2")
for q in (1.5, 2, 3, 5, 7):
    k0 = int(round(np.log(q) / du)); rho = [1 - (P(g, k0) / p0) / Rg[k0] for g, p0 in zip(zs, P0)]
    print(f"   {q:4.1f}   {np.log10(Rg[k0]):+7.3f}          " + "   ".join(f"{x:.3f}" for x in rho) + "                    " + "   ".join(f"{x:.3f}" for x in meas[q]))
