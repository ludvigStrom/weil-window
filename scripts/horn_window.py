"""The horn's axis cut at height T, against the full wave (notes/picture.md round 5: comparing the two cutoffs).
On the axis the wave of energy eps is F(y) = A y^s + B y^(1-s) + N(y), F(1/y) = F(y). A window [1/T, T] on the axis is the Weil
window [-a, a] with T = e^a = sqrt(c). Its transform in scaling frequencies is
    W_T(t) = int_{1/T}^{T} F(y) y^(it) dy/y = 2 int_1^T F(y) cos(t ln y) dy/y.
Exact decomposition (derived in the notes):   W_T(t) = xi(1/2 + i(t+eps)) xi(1/2 + i(t-eps)) + Wall_T(t) - Tail_T(t),
    Wall_T(t) = A T^(s+it)/(s+it) + A T^(s-it)/(s-it) + B T^(1-s+it)/(1-s+it) + B T^(1-s-it)/(1-s-it)   (the wave's two plane parts at the cut)
    Tail_T(t) = 2 int_T^inf N(y) cos(t ln y) dy/y.
Printed: a check of the identity, the sizes of the three pieces, and the real zeros of W_T against those of the wall term alone.
Usage: python3 horn_window.py eps c tmax"""
import sys, time, math
from mpmath import mp, mpf, mpc, besselk, pi, sqrt, log, exp, gamma, zeta, zetazero
from mpmath.calculus.quadrature import GaussLegendre
mp.dps = 50; eps = mpf(sys.argv[1]); c = mpf(sys.argv[2]); tmax = float(sys.argv[3]); T = sqrt(c); s = mpc(0.5, eps); t0 = time.time()
xi = lambda z: pi ** (-z / 2) * gamma(z / 2) * zeta(z); A = xi(1 + 2j * eps); B = xi(2j * eps)
tau = [None] + [sum(mp.cos(eps * log(mpf(a) * a / n)) for a in range(1, n + 1) if n % a == 0) for n in range(1, 41)]
def Nw(y):
    tot = mpf(0)
    for n in range(1, 41):
        x = 2 * pi * n * y
        if x > 130: break
        tot += tau[n] * besselk(mpc(0, eps), x).real
    return 4 * sqrt(y) * tot
nodes = GaussLegendre(mp).calc_nodes(5, mp.prec)
def grid(u0, u1, P):
    g = []
    for k in range(P):
        a, b = u0 + (u1 - u0) * k / P, u0 + (u1 - u0) * (k + 1) / P
        for x, w in nodes: u = (a + b) / 2 + (b - a) / 2 * x; g.append((u, w * (b - a) / 2, Nw(exp(u))))
    return g
gin = grid(mpf(0), log(T), 8); gout = grid(log(T), log(mpf(14)), 10)
def pieces(t):
    t = mpf(t); win = sum(wt * Nv * 2 * mp.cos(t * u) for u, wt, Nv in gin); tail = sum(wt * Nv * 2 * mp.cos(t * u) for u, wt, Nv in gout)
    cT = sum(A * (T ** (s + sg * 1j * t) - 1) / (s + sg * 1j * t) + B * (T ** (1 - s + sg * 1j * t) - 1) / (1 - s + sg * 1j * t) for sg in (1, -1))
    wall = sum(A * T ** (s + sg * 1j * t) / (s + sg * 1j * t) + B * T ** (1 - s + sg * 1j * t) / (1 - s + sg * 1j * t) for sg in (1, -1))
    return (win + cT).real, (xi(mpc(0.5, t + eps)) * xi(mpc(0.5, t - eps))).real, wall.real, tail
print(f"eps = {eps}, window c = {c} (cut at height T = {mp.nstr(T, 5)}), {len(gin) + len(gout)} heights  ({time.time() - t0:.0f} s)")
print("   t | window W_T | zeta x zeta | wall term | tail | W_T - (zz + wall - tail)")
for t in (0, 3, 6, 10, 14.1347, 20, 30):
    W, zz, wall, tail = pieces(t); print(f"   {t:7.4f} | {mp.nstr(W, 8):>14s} | {mp.nstr(zz, 4):>10s} | {mp.nstr(wall, 8):>14s} | {mp.nstr(tail, 3):>10s} | {mp.nstr(W - (zz + wall - tail), 3)}")
def zeros(f, n=int(tmax * 20)):
    out = []; prev = f(0)
    for i in range(1, n + 1):
        t = tmax * i / n; v = f(t)
        if prev * v < 0: out.append(t - (tmax / n) * v / (v - prev))
        prev = v
    return out
zW = zeros(lambda t: pieces(t)[0]); zwall = zeros(lambda t: pieces(t)[2])
print(f"   real zeros of the window's transform up to t = {tmax}: {len(zW)};  of the wall term alone: {len(zwall)};  window length / pi times tmax = {float(2 * log(T)) / math.pi * tmax:.1f}")
print("   window: " + " ".join(f"{float(z):.2f}" for z in zW[:14])); print("   wall:   " + " ".join(f"{float(z):.2f}" for z in zwall[:14]))
print("   zeta zeros +- eps would be: " + " ".join(f"{float(zetazero(k).imag) + sg * float(eps):.2f}" for k in (1, 2, 3) for sg in (-1, 1)))
print(f"   ({time.time() - t0:.0f} s)")
