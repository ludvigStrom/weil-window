"""The multiplier K_n = (tick_n - zero_n) / (vacuum level)  (notes/questions.md, 'The K multiplier').
Pairing Q v0 = lambda0 v0 (Q = 2 sum over zeros |e_g><e_g|) with the vacuum-minus-one-zero wave gives the identity
   K_n = (t_n / v0^'(t_n)^2) [ P_n(0) + 2 sum_{g != g_n} w_g / (g^2 - t_n^2) ],   w_g = v0^(g)^2 / lambda0,  sum 2 w_g = 1,
with v0^ the cosine transform of the normalised vacuum and P_n(0) = (1/g_n) int_0^inf R(y) sin(g_n y) dy.
First part in the limit of a large window, from Riemann's Phi alone:  K_n^Phi = g_n P_n(0) / Xi~'(g_n)^2,
Xi~ = transform of Phi / |Phi|.  Compared with K_n measured on windows.   Usage: python3 k_multiplier.py"""
import numpy as np
du = 2.5e-4; u = np.arange(-4, 4 + du / 2, du); n = np.arange(1, 7)[:, None]
def Phi(x):
    x = np.abs(x); return 2 * np.sum((2 * np.pi ** 2 * n ** 4 * np.exp(4.5 * x) - 3 * np.pi * n ** 2 * np.exp(2.5 * x)) * np.exp(-np.pi * n ** 2 * np.exp(2 * x)), axis=0)
F = Phi(u); F /= np.sqrt(np.sum(F * F) * du); M = int(round(2.6 / du)); Rg = np.array([np.sum(F[:len(F) - k] * F[k:]) * du for k in range(M + 1)]); yg = np.arange(M + 1) * du
zs = {1: 14.134725141734695, 2: 21.022039638771556, 3: 25.010857580145689, 4: 30.424876125859512}
meas = {14.2: {1: 2.07e-61 / 2.817e-65}, 18.2: {1: 2.98e-81 / 3.571e-85}, 22.2: {1: 6.88e-100 / 7.555e-104}, 28.2: {1: 9923.91, 2: 2.88514e7, 3: 3.76445e9, 4: 5.81843e12}}
print("zero    g^2 P(0)    Xi~'(g)        K^Phi (first part, window -> infinity)     K measured at c = 28.2     ratio")
for k, g in zs.items():
    P0 = np.sum(Rg * np.sin(g * yg)) * du / g; dX = -np.sum(u * F * np.sin(g * u)) * du; K = g * P0 / dX ** 2
    print(f"  {k}     {g * g * P0:.4f}     {dX:+.4e}      {K:.5e}                              {meas[28.2][k]:.5e}          {meas[28.2][k] / K:.4f}")
g = zs[1]; P0 = np.sum(Rg * np.sin(g * yg)) * du / g; dX = -np.sum(u * F * np.sin(g * u)) * du; K = g * P0 / dX ** 2
print("first zero, K measured / K^Phi by window: " + ", ".join(f"c = {c}: {m[1] / K:.4f}" for c, m in meas.items()))
