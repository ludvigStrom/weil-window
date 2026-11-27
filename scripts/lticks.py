"""The recipe for a Dirichlet L-function: lowest eigenvalues of both halves of the window form, and the real zeros of the
transform of the lowest even state against the zeros of L(1/2 + i t, chi).   Usage: python3 lticks.py chi4 N c SMAX   (env DPS)"""
import sys, math, os
import mpmath
from mpmath import mp
from flint import arb
import hp_block as hb, lfunc_block as lf
mp.dps = int(os.environ.get("DPS", 60)); name = sys.argv[1]; N = int(sys.argv[2]); c = float(sys.argv[3]); SMAX = float(sys.argv[4])
F = lf.FormChi(c, name); L = F.L; tomp = lambda x: mp.mpf(x.mid().str(hb.DIG, radius=False)); Lm = tomp(L); om = [2 * mp.pi * m / Lm for m in range(N + 1)]
nodes = hb.quad_nodes(L, [arb(q).log() for q, p in F.pp], 2 * math.pi * N / float(L.mid())); S_, C2_ = hb.transforms(nodes, [F.H(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, F.H(hb.ZERO), F.H(L), L, N)
def eig(Q, n):
    M = mp.matrix(n, n)
    for i in range(n):
        for j in range(n): M[i, j] = (tomp(Q[i, j]) + tomp(Q[j, i])) / 2
    return mp.eigsy(M)
Ee, Ue = eig(Qe, N + 1); Eo, Uo = eig(Qo, N); se = sorted(Ee); so = sorted(Eo)
print(f"{name} (conductor {F.Nc}, parity {F.a}), c = {c}: prime powers with their signs: " + ", ".join(f"{'+' if F.chi(q) > 0 else '-'}{q}" for q, p in F.pp))
print("   lowest even eigenvalues: " + ", ".join(mp.nstr(x, 4) for x in se[:5])); print("   lowest odd eigenvalues:  " + ", ".join(mp.nstr(x, 4) for x in so[:5]))
k = min(range(N + 1), key=lambda i: Ee[i]); v = [Ue[i, k] for i in range(N + 1)]
def xi_hat(s): return sum(v[m] * (1 / mp.sqrt(Lm) if m == 0 else mp.sqrt(2 / Lm)) * (mp.sincpi((s - om[m]) * Lm / (2 * mp.pi)) + mp.sincpi((s + om[m]) * Lm / (2 * mp.pi))) for m in range(N + 1))
ss = [mp.mpf(i) / 20 for i in range(1, int(SMAX * 20) + 1)]; vals = [xi_hat(s) for s in ss]; ticks = [mp.findroot(xi_hat, (ss[i], ss[i + 1]), solver="anderson") for i in range(len(ss) - 1) if vals[i] * vals[i + 1] < 0]
tab = [0] * F.Nc
for r in range(F.Nc): tab[r] = F.chi(r)
Lf = lambda t: mp.dirichlet(mp.mpf(1) / 2 + 1j * t, tab); zs = []
for t in ticks:
    try: z = mp.findroot(Lf, mp.mpc(t, 0), tol=mp.mpf(10) ** (-mp.dps + 8)); zs.append(z)
    except Exception: zs.append(None)
print(f"   {len(ticks)} real ticks below {SMAX:g}: " + ", ".join(mp.nstr(t, 10) for t in ticks))
print("   nearest zero of L(1/2 + it): " + ", ".join("none" if z is None else mp.nstr(mp.re(z), 10) for z in zs))
print("   tick - zero: " + ", ".join("-" if z is None else mp.nstr(t - mp.re(z), 2) for t, z in zip(ticks, zs)))
