"""The escaping wave of the modular horn, built from divisor sums alone, and its absorption lines (notes/picture.md, round 5).
State of energy eps (s = 1/2 + i eps), on the axis x = 0 of the horn:
    F(y) = A y^s + B y^(1-s) + N(y),   N(y) = 4 sqrt(y) sum_n tau(n) K_{i eps}(2 pi n y),   tau(n) = sum_{ab = n} (a/b)^(i eps).
Nothing about zeta is put in. A and B (the incoming and outgoing amplitudes) are found by asking that the wave fits the horn,
F(1/y) = F(y), at two heights, and checked at a third. Then the wave is written in scaling frequencies t:
    M(1/2 + it) = int_0^inf N(y) y^(it) dy/y,
done by folding at y = 1 (the wave's own symmetry), so only y >= 1 is needed. Prediction (Hecke, Zagier):
    M(w) = xi(w + i eps) xi(w - i eps),  xi(s) = pi^(-s/2) Gamma(s/2) zeta(s):  zeros at t = gamma +- eps for every zeta zero gamma;
    B / A = xi(2 i eps) / xi(1 + 2 i eps)  (the horn's reflection coefficient).
Usage: python3 horn_wave.py eps tmax"""
import sys, json, time
from mpmath import mp, mpf, mpc, besselk, pi, sqrt, log, exp, gamma, zeta, matrix, lu_solve, zetazero
from mpmath.calculus.quadrature import GaussLegendre
mp.dps = 80; eps = mpf(sys.argv[1]); tmax = float(sys.argv[2]); t0 = time.time(); s = mpc(0.5, eps); NMAX = 60
tau = [None] + [sum(mp.cos(eps * log(mpf(a) * a / n)) for a in range(1, n + 1) if n % a == 0) for n in range(1, NMAX + 1)]
def Nw(y):
    tot = mpf(0)
    for n in range(1, NMAX + 1):
        x = 2 * pi * n * y
        if x > 200: break
        tot += tau[n] * besselk(mpc(0, eps), x).real
    return 4 * sqrt(y) * tot
def solveAB(y1, y2):
    M = matrix(2, 2); r = matrix(2, 1)
    for i, y in enumerate((y1, y2)):
        M[i, 0] = y ** s - y ** (-s); M[i, 1] = y ** (1 - s) - y ** (s - 1); r[i] = Nw(1 / y) - Nw(y)
    return lu_solve(M, r)
ab = solveAB(mpf("1.25"), mpf("1.6")); A, B = ab[0], ab[1]; y3 = mpf("1.4")
res = A * (y3 ** s - y3 ** (-s)) + B * (y3 ** (1 - s) - y3 ** (s - 1)) - (Nw(1 / y3) - Nw(y3))
xi = lambda z: pi ** (-z / 2) * gamma(z / 2) * zeta(z)
print(f"eps = {eps}: wave fitted to the horn at two heights; misfit at a third height {mp.nstr(abs(res), 3)}  ({time.time() - t0:.0f} s)")
print(f"  reflection coefficient B/A from the wave: {mp.nstr(B / A, 12)}"); print(f"  xi(2 i eps)/xi(1 + 2 i eps):               {mp.nstr(xi(2j * eps) / xi(1 + 2j * eps), 12)}   (modulus {mp.nstr(abs(B / A), 12)})")
# fold at y = 1: M(w) = int_1^inf N(y) [y^(w-1/2) + y^(1/2-w)] dy/y + A/(w-1/2-s) + B/(w-3/2+s) - A/(w-1/2+s) - B/(w+1/2-s)
nodes = GaussLegendre(mp).calc_nodes(5, mp.prec); U = log(mpf(16)); P = 24; grid = []
for k in range(P):
    a, b = U * k / P, U * (k + 1) / P
    for x, w in nodes: u = (a + b) / 2 + (b - a) / 2 * x; grid.append((u, w * (b - a) / 2, Nw(exp(u))))
print(f"  wave tabulated at {len(grid)} heights  ({time.time() - t0:.0f} s)")
def M(t):
    w = mpc(0.5, t); I = sum(wt * Nv * 2 * mp.cos(t * u) for u, wt, Nv in grid)        # y^(it) + y^(-it) = 2 cos(t ln y)
    return I + A / (w - 0.5 - s) + B / (w - 1.5 + s) - A / (w - 0.5 + s) - B / (w + 0.5 - s)
out = []; nt = int(tmax * 5) + 1
for i in range(nt):
    t = mpf(i) / 5; m = M(t); norm = abs(pi ** (-mpc(0.5, t)) * gamma(mpc(0.25, (t + eps) / 2)) * gamma(mpc(0.25, (t - eps) / 2)))
    exact = xi(mpc(0.5, t + eps)) * xi(mpc(0.5, t - eps))
    out.append(dict(t=float(t), wave=float(m.real / norm), imag=float(abs(m.imag) / norm), formula=float(exact.real / norm)))
err = max(abs(o["wave"] - o["formula"]) for o in out); print(f"  largest difference between the wave's transform and xi xi (after dividing out the Gamma factors): {err:.2e};  largest imaginary part {max(o['imag'] for o in out):.2e}")
# zeros of the wave's transform by sign changes, against gamma +- eps
zs = []
for a, b in zip(out[:-1], out[1:]):
    if a["wave"] * b["wave"] < 0: zs.append(a["t"] - a["wave"] * (b["t"] - a["t"]) / (b["wave"] - a["wave"]))
gam = [float(zetazero(k).imag) for k in range(1, 12)]; pred = sorted(x for g in gam for x in (g - float(eps), g + float(eps)) if 0 < x < tmax)
print("  dark lines found in the wave:      " + " ".join(f"{z:.2f}" for z in zs)); print("  zeta zeros plus and minus eps:     " + " ".join(f"{z:.2f}" for z in pred))
json.dump(dict(eps=float(eps), A=[float(A.real), float(A.imag)], B=[float(B.real), float(B.imag)], curve=out, zeros=zs, predicted=pred), open(f"horn_wave_{float(eps)}.json", "w"))
print(f"  ({time.time() - t0:.0f} s)")
