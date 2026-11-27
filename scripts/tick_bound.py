"""A candidate bound on how far a tick is from its limit, built from the window's own data (notes/questions.md, 'A bound to aim for').
For the lowest even state v (unit norm, Fourier coefficients v_m on frequencies w_m), its cosine transform v^, its n-th tick t_n:
      B_n(c) = t_n * P_n * eps / v^'(t_n)^2 ,    P_n = sum_m v_m^2 / (t_n^2 - w_m^2) ,   eps = lowest eigenvalue.
Everything in B_n comes from the form on the window c; no zero of zeta or L enters. Pairing the eigen-equation with v divided by
(s^2 - t_n^2) gives   t_n - gamma_n = B_n (1 + r_n),   r_n = 2 sum_{gamma != gamma_n} w_gamma / ((gamma^2 - t_n^2) P_n),
so the conjecture to test is that  rho_n = (t_n - gamma_n) / B_n  stays between fixed constants (near 1) at every window.
Usage: python3 tick_bound.py N c1,c2,...  [chi]    env DPS (default 130).  Appends to tick_bound.json"""
import sys, math, os, json, time
import mpmath
from mpmath import mp
from flint import arb
import hp_block as hb
mp.dps = int(os.environ.get("DPS", 130)); N = int(sys.argv[1]); CS = [float(x) for x in sys.argv[2].split(",")]; CHI = sys.argv[3] if len(sys.argv) > 3 else None; NT = int(os.environ.get("NT", 4))
if CHI: import lfunc_block as lf
tomp = lambda x: mp.mpf(x.mid().str(hb.DIG, radius=False)); fn = "tick_bound.json"; zc = "tick_bound_zeros.json"; Z = json.load(open(zc)) if os.path.exists(zc) else {}
def true_zeros(guess):
    key = (CHI or "zeta") + f"|{mp.dps}"
    if key not in Z:
        if CHI:
            Nc = lf.CHARS[CHI][0]; tab = [lf.CHARS[CHI][2].get(r, 0) for r in range(Nc)]; Lf = lambda t: mp.dirichlet(mp.mpf(1) / 2 + 1j * t, tab)
            Z[key] = [str(mp.re(mp.findroot(Lf, mp.mpc(g, 0), tol=mp.mpf(10) ** (-mp.dps + 10)))) for g in guess]
        else: Z[key] = [str(mp.im(mp.zetazero(n + 1))) for n in range(len(guess))]
        json.dump(Z, open(zc, "w"))
    return [mp.mpf(x) for x in Z[key]]
for c in CS:
    T0 = time.time(); F = lf.FormChi(c, CHI) if CHI else hb.Form(c); L = F.L; Lm = tomp(L); n1 = N + 1; om = [2 * mp.pi * m / Lm for m in range(n1)]; nrm = [1 / mp.sqrt(Lm)] + [mp.sqrt(2 / Lm)] * N
    nodes = hb.quad_nodes(L, [arb(q).log() for q, pp in F.pp], 2 * math.pi * N / float(L.mid())); S_, C2_ = hb.transforms(nodes, [F.H(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, F.H(hb.ZERO), F.H(L), L, N)
    Q = mp.matrix(n1, n1)
    for i in range(n1):
        for j in range(n1): Q[i, j] = (tomp(Qe[i, j]) + tomp(Qe[j, i])) / 2
    Qi = mp.inverse(Q); v = mp.matrix([mp.mpf(1) / (i + 1) for i in range(n1)])
    for it in range(6): v = Qi * v; v = v / mp.norm(v)
    eps = 1 / mp.norm(Qi * v); hat = lambda s: sum(v[m] * nrm[m] * (mp.sincpi((s - om[m]) * Lm / (2 * mp.pi)) + mp.sincpi((s + om[m]) * Lm / (2 * mp.pi))) for m in range(n1)) * Lm / 2
    ss = [mp.mpf(i) / 10 for i in range(10, 400)]; vals = []; ticks = []
    for i, s in enumerate(ss):
        vals.append(hat(s))
        if i and vals[-1] * vals[-2] < 0: ticks.append(mp.findroot(hat, (ss[i - 1], s), solver="anderson"))
        if len(ticks) == NT: break
    gam = true_zeros([float(t) for t in ticks]); row = dict(c=c, N=N, chi=CHI or "zeta", eps=float(mp.log10(eps)), n=[])
    for n, (t, g) in enumerate(zip(ticks, gam)):
        d1 = mp.diff(hat, t); P = sum(v[m] ** 2 / (t * t - om[m] ** 2) for m in range(n1)); B = t * P * eps / d1 ** 2; rho = (t - g) / B
        row["n"].append(dict(tick=float(t), off=float(mp.log10(abs(t - g))), B=float(mp.log10(abs(B))), rho=float(rho), t2P=float(t * t * P)))
    print(f"{CHI or 'zeta'} c = {c:6.2f}: eps = 1e{float(mp.log10(eps)):.2f}; " + "; ".join(f"tick {i + 1}: off 1e{r['off']:.2f}, bound 1e{r['B']:.2f}, rho {r['rho']:.4f}" for i, r in enumerate(row["n"])) + f"   ({time.time() - T0:.0f} s)", flush=True)
    R = json.load(open(fn)) if os.path.exists(fn) else []; R.append(row); json.dump(R, open(fn, "w"))
