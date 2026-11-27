"""A prime entering the window (notes/questions.md, 'the kink'). For windows c just above a prime p:
the share of the vacuum's level carried by that prime's mirror,  E = w <v0|A_ln p|v0> / lambda_0  (minus d ln lambda_0 / d ln w, exact),
against h = ln(c / p), the width of the strip the mirror connects. If the vacuum had a finite value at the edge, E would grow like h;
if it vanishes like 1 / sqrt(ln(1 / (kappa t))) (string.md), E grows like h / ln(1 / (kappa h)).
Usage: python3 prime_entry.py N c p"""
import sys, math, os, json, time
import mpmath
from mpmath import mp
from flint import arb
import hp_block as hb
mp.dps = int(os.environ.get("DPS", 170)); N = int(sys.argv[1]); c = float(sys.argv[2]); T0 = time.time()
F = hb.Form(c); L = F.L; tomp = lambda x: mp.mpf(x.mid().str(hb.DIG, radius=False)); Lm = tomp(L); n1 = N + 1
om = [2 * mp.pi * m / Lm for m in range(2 * N + 2)]; nrm = [1 / mp.sqrt(Lm)] + [mp.sqrt(2 / Lm)] * N
p = int(sys.argv[3])
nodes = hb.quad_nodes(L, [arb(q).log() for q, pp in F.pp], 2 * math.pi * N / float(L.mid()))
S_, C2_ = hb.transforms(nodes, [F.H(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, F.H(hb.ZERO), F.H(L), L, N)
Q = mp.matrix(n1, n1)
for i in range(n1):
    for j in range(n1): Q[i, j] = (tomp(Qe[i, j]) + tomp(Qe[j, i])) / 2
print(f"c = {c}, N = {N}: form built ({time.time() - T0:.0f} s)", flush=True); Qi = mp.inverse(Q); v = mp.matrix([mp.mpf(1) / (i + 1) for i in range(n1)])
for it in range(6): v = Qi * v; v = v / mp.norm(v)
lam0 = (v.T * Q * v)[0]; w2 = Qi * v; print(f"   vacuum eigenvalue {mp.nstr(lam0, 4)} (check: 1 / |Q^-1 v| = {mp.nstr(1 / mp.norm(w2), 4)}); inverse done ({time.time() - T0:.0f} s)", flush=True)
def Aops(a):
    """A_a and dA_a/da on the even basis, closed form."""
    lo, hi = -Lm / 2, Lm / 2 - a; G = mp.matrix(n1, n1); D = mp.matrix(n1, n1)
    sl = [mp.sin(o * lo) for o in om]; cl = [mp.cos(o * lo) for o in om]; sh = [mp.sin(o * hi) for o in om]; ch = [mp.cos(o * hi) for o in om]; sa = [mp.sin(om[n] * a) for n in range(n1)]; ca = [mp.cos(om[n] * a) for n in range(n1)]; cL = [mp.cos(om[n] * Lm / 2) for n in range(n1)]
    def tab(k):                                    # sin, cos of om_k * x at x = hi, lo for signed index k
        s = 1 if k >= 0 else -1; k = abs(k); return s * sh[k], ch[k], s * sl[k], cl[k]
    for m in range(n1):
        for n in range(n1):
            tot_c = tot_s = mp.mpf(0)
            for k, sgn in ((m - n, -1), (m + n, 1)):                # integrand cos / sin (om_k u + sgn om_n a)
                shk, chk, slk, clk = tab(k); sp, cp = sgn * sa[n], ca[n]
                if k == 0: Ic = (hi - lo) * cp; Is = (hi - lo) * sp
                else:
                    kap = om[abs(k)] * (1 if k > 0 else -1)
                    Ic = ((shk * cp + chk * sp) - (slk * cp + clk * sp)) / kap                     # int cos(kap u + phi)
                    Is = -((chk * cp - shk * sp) - (clk * cp - slk * sp)) / kap                    # int sin(kap u + phi)
                tot_c += Ic; tot_s += sgn * Is
            G[m, n] = nrm[m] * nrm[n] * tot_c / 2; D[m, n] = nrm[m] * nrm[n] * (-om[n] * tot_s / 2 - ch[m] * cL[n])
    return G + G.T, D + D.T
a = mp.log(p); A, dA = Aops(a); w = mp.log(p) / mp.sqrt(p); E = w * (v.T * A * v)[0] / lam0; h = Lm - a; kap = 2 * mp.pi * mp.exp(mp.euler)
print(f"RESULT c = {c}: h = {float(h):.5f}; level {mp.nstr(lam0, 5)}; share of the level from prime {p}: E = {mp.nstr(E, 6)}; E / h = {mp.nstr(E / h, 6)}; E ln(1/(kappa h)) / h = {mp.nstr(E * mp.log(1 / (kap * h)) / h, 6)}", flush=True)
