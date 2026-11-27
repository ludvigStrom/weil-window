"""A fingerprint of passivity (notes/questions.md, 'Passive system'). For a lossless passive system the resonances and
anti-resonances alternate (Foster; Hermite-Biehler for de Branges spaces). Two natural sets of frequencies on a window:
the real zeros of the cosine transform of the lowest even state (the ticks, near the zeros of zeta) and the real zeros of
the sine transform of the lowest odd state. Do they interlace?   Usage: python3 interlace.py N c SMAX   (env DPS)"""
import sys, math, os
import mpmath
from mpmath import mp
from flint import arb
import hp_block as hb
mp.dps = int(os.environ.get("DPS", 60)); N = int(sys.argv[1]); c = float(sys.argv[2]); SMAX = float(sys.argv[3]); CHI = os.environ.get("CHI")
if CHI:
    import lfunc_block as lf; F = lf.FormChi(c, CHI)
else: F = hb.Form(c)
L = F.L; tomp = lambda x: mp.mpf(x.mid().str(hb.DIG, radius=False)); Lm = tomp(L); om = [2 * mp.pi * m / Lm for m in range(N + 1)]
nodes = hb.quad_nodes(L, [arb(q).log() for q, p in F.pp], 2 * math.pi * N / float(L.mid())); S_, C2_ = hb.transforms(nodes, [F.H(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, F.H(hb.ZERO), F.H(L), L, N)
def low(Q, n):
    M = mp.matrix(n, n)
    for i in range(n):
        for j in range(n): M[i, j] = (tomp(Q[i, j]) + tomp(Q[j, i])) / 2
    E, U = mp.eigsy(M); k = min(range(n), key=lambda i: E[i]); return E[k], [U[i, k] for i in range(n)]
ee, ve = low(Qe, N + 1); eo, vo = low(Qo, N); sc = lambda t: mp.sincpi(t * Lm / (2 * mp.pi))
fe = lambda s: sum(ve[m] * (1 / mp.sqrt(Lm) if m == 0 else mp.sqrt(2 / Lm)) * (sc(s - om[m]) + sc(s + om[m])) for m in range(N + 1))
fo = lambda s: sum(vo[m - 1] * mp.sqrt(2 / Lm) * (sc(s - om[m]) - sc(s + om[m])) for m in range(1, N + 1))
ss = [mp.mpf(i) / 20 for i in range(1, int(SMAX * 20) + 1)]
def zeros(f):
    v = [f(s) for s in ss]; return [mp.findroot(f, (ss[i], ss[i + 1]), solver="anderson") for i in range(len(ss) - 1) if v[i] * v[i + 1] < 0]
te, to = zeros(fe), zeros(fo); allz = sorted([(t, "E") for t in te] + [(t, "O") for t in to]); pat = "".join(k for t, k in allz)
print(f"{CHI or 'zeta'}, c = {c}: lowest eigenvalue even {mp.nstr(ee, 3)}, odd {mp.nstr(eo, 3)}")
print("   even vacuum, zeros of the cosine transform: " + ", ".join(mp.nstr(t, 8) for t in te)); print("   odd vacuum, zeros of the sine transform:    " + ", ".join(mp.nstr(t, 8) for t in to))
print(f"   order along the axis: {pat}   -> " + ("strictly alternating" if all(pat[i] != pat[i + 1] for i in range(len(pat) - 1)) else "NOT alternating"))
if te and to: print("   smallest distance between an even and an odd zero: " + mp.nstr(min(abs(a - b) for a in te for b in to), 3))
