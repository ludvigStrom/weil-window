"""Connes's recipe, and what happens to it when the primes are moved (notes/questions.md, 'The clock's ticks').
Lowest even eigenvector xi of the form on the window (in mpmath, DPS digits); its transform
xi^(s) = int xi(u) cos(s u) du is an even entire function; its real zeros are the 'ticks'
(the spectrum of the scaling operator with the boundary rule built from xi; CCM arXiv:2511.22755).
Moves: zoom p -> p^(1 + eps) (ZOOM), or one prime p0 -> p0 (1 + delta).
Usage: python3 clock_ticks.py N c SMAX [zoom eps | move p0 delta]      env DPS (default 60)"""
import sys, math, os, json
import mpmath
from mpmath import mp
from flint import arb
import hp_block as hb
mp.dps = int(os.environ.get("DPS", 60)); N = int(sys.argv[1]); c = float(sys.argv[2]); SMAX = float(sys.argv[3]); mode = sys.argv[4] if len(sys.argv) > 4 else "true"; par = sys.argv[5:]
F = hb.Form(c); L = F.L; primes = sorted({pp for q, pp in F.pp}); ALL = tuple(q for q, pp in F.pp); tomp = lambda x: mp.mpf(x.mid().str(hb.DIG, radius=False))
if mode == "zoom": e = arb(par[0]); xs = [arb(p) ** (1 + e) for p in primes]
elif mode == "move": xs = [arb(p) * (1 + arb(par[1])) if p == int(par[0]) else arb(p) for p in primes]
else: xs = [arb(p) for p in primes]
def even_matrix(xs):
    terms = [(x.log(), x ** k) for x in xs for k in range(1, 40) if float((x ** k).mid()) < c]
    Hf = lambda U: F.H(U, exclude=ALL) + sum((lg * F.ramp(U, q) for lg, q in terms), hb.ZERO)
    nodes = hb.quad_nodes(L, [q.log() for lg, q in terms], 2 * math.pi * N / float(L.mid()))
    S_, C2_ = hb.transforms(nodes, [Hf(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, Hf(hb.ZERO), Hf(L), L, N)
    M = mp.matrix(N + 1, N + 1)
    for i in range(N + 1):
        for j in range(N + 1): M[i, j] = (tomp(Qe[i, j]) + tomp(Qe[j, i])) / 2
    return M
E, U = mp.eigsy(even_matrix(xs)); k = min(range(N + 1), key=lambda i: E[i]); FOLLOW = os.environ.get("FOLLOW") and mode != "true"
if FOLLOW:                                         # follow the true vacuum instead of the lowest state
    E0, U0 = mp.eigsy(even_matrix([arb(p) for p in primes])); k0 = min(range(N + 1), key=lambda i: E0[i]); v0 = [U0[i, k0] for i in range(N + 1)]
    ov = [abs(sum(v0[i] * U[i, j] for i in range(N + 1))) for j in range(N + 1)]; k = max(range(N + 1), key=lambda j: ov[j]); rank = sorted(E).index(E[k]) + 1
    print(f"   following the true vacuum: best match is eigenstate number {rank} from the bottom, eigenvalue {mp.nstr(E[k], 4)}, overlap {mp.nstr(ov[k], 6)}; next best overlap {mp.nstr(sorted(ov)[-2], 3)}")
v = [U[i, k] for i in range(N + 1)]; Lm = tomp(L); om = [2 * mp.pi * m / Lm for m in range(N + 1)]
srt = sorted(E); print(f"c = {c}, {mode} {' '.join(par)}: lowest even eigenvalues " + ", ".join(mp.nstr(x, 4) for x in srt[:4]))
def xi_hat(s):
    tot = mp.mpf(0)
    for m in range(N + 1):
        nrm = 1 / mp.sqrt(Lm) if m == 0 else mp.sqrt(2 / Lm)
        tot += v[m] * nrm * (mp.sincpi((s - om[m]) * Lm / (2 * mp.pi)) + mp.sincpi((s + om[m]) * Lm / (2 * mp.pi))) * Lm / 2
    return tot
ss = [mp.mpf(i) / 20 for i in range(1, int(SMAX * 20) + 1)]; vals = [xi_hat(s) for s in ss]; ticks = []
for i in range(len(ss) - 1):
    if vals[i] * vals[i + 1] < 0: ticks.append(mp.findroot(xi_hat, (ss[i], ss[i + 1]), solver="anderson"))
zz = []
n = 1
while True:
    z = mp.im(mp.zetazero(n))
    if z > SMAX: break
    zz.append(z); n += 1
scale = 1 + mp.mpf(par[0]) if mode == "zoom" else mp.mpf(1)
print(f"   {len(ticks)} real ticks below {SMAX:g}; zeta has {len(zz)} zeros there")
print("   ticks:      " + ", ".join(mp.nstr(t, 10) for t in ticks))
near = lambda t, ref: min(ref, key=lambda z: abs(z - t))
print("   tick - nearest zero:           " + ", ".join(mp.nstr(t - near(t, zz), 2) for t in ticks))
if mode == "zoom": print("   tick - nearest zero/(1+eps):   " + ", ".join(mp.nstr(t - near(t, [z / scale for z in zz]), 2) for t in ticks))
json.dump(dict(c=c, N=N, mode=mode, par=par, ticks=[float(t) for t in ticks], eig=[float(x) for x in srt[:4]]), open(f"clock_ticks_{c:g}_{mode}{'_' + '_'.join(par) if par else ''}{'_follow' if FOLLOW else ''}.json", "w"))
