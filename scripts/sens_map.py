"""The map: how much does each zero care where each prime is? (notes/questions.md, 'The map')
Ticks = real zeros of the transform of the lowest even state (clock_ticks.py). Each prime p is moved to p (1 + d) with d so
small that the vacuum stays the lowest state and responds linearly (d = DELTA, default 1e-70; one-sided difference; env PRIMES=2,3 does only those and merges into the json).
Cell = d(tick)/d(delta). The 'tolerance' of a zero to a prime is then (precision wanted) / cell.
Usage: python3 sens_map.py N c SMAX     env DPS (default 90), DELTA"""
import sys, math, os, json
import mpmath
from mpmath import mp
from flint import arb
import hp_block as hb
mp.dps = int(os.environ.get("DPS", 150)); DELTA = os.environ.get("DELTA", "1e-70"); ONLY = [int(x) for x in os.environ["PRIMES"].split(",")] if os.environ.get("PRIMES") else None; N = int(sys.argv[1]); c = float(sys.argv[2]); SMAX = float(sys.argv[3])
F = hb.Form(c); L = F.L; primes = sorted({pp for q, pp in F.pp}); ALL = tuple(q for q, pp in F.pp); tomp = lambda x: mp.mpf(x.mid().str(hb.DIG, radius=False)); Lm = tomp(L); om = [2 * mp.pi * m / Lm for m in range(N + 1)]
def vacuum(xs):
    terms = [(x.log(), x ** k) for x in xs for k in range(1, 40) if float((x ** k).mid()) < c]
    Hf = lambda U: F.H(U, exclude=ALL) + sum((lg * F.ramp(U, q) for lg, q in terms), hb.ZERO)
    nodes = hb.quad_nodes(L, [q.log() for lg, q in terms], 2 * math.pi * N / float(L.mid()))
    S_, C2_ = hb.transforms(nodes, [Hf(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, Hf(hb.ZERO), Hf(L), L, N); M = mp.matrix(N + 1, N + 1)
    for i in range(N + 1):
        for j in range(N + 1): M[i, j] = (tomp(Qe[i, j]) + tomp(Qe[j, i])) / 2
    E, U = mp.eigsy(M); k = min(range(N + 1), key=lambda i: E[i]); return [U[i, k] for i in range(N + 1)], sorted(E)[:3]
def xi_hat(v):
    def f(s):
        tot = mp.mpf(0)
        for m in range(N + 1):
            nrm = 1 / mp.sqrt(Lm) if m == 0 else mp.sqrt(2 / Lm); tot += v[m] * nrm * (mp.sincpi((s - om[m]) * Lm / (2 * mp.pi)) + mp.sincpi((s + om[m]) * Lm / (2 * mp.pi)))
        return tot
    return f
v0, e0 = vacuum([arb(p) for p in primes]); f0 = xi_hat(v0); ss = [mp.mpf(i) / 10 for i in range(1, int(SMAX * 10) + 1)]; vals = [f0(s) for s in ss]
t0 = [mp.findroot(f0, (ss[i], ss[i + 1]), solver="anderson") for i in range(len(ss) - 1) if vals[i] * vals[i + 1] < 0]
print(f"c = {c}: primes {primes}; lowest even eigenvalues " + ", ".join(mp.nstr(x, 3) for x in e0) + f"; {len(t0)} ticks below {SMAX:g}; delta = {DELTA}", flush=True)
fn = f"sens_map_{c:g}.json"; out = json.load(open(fn)) if os.path.exists(fn) and ONLY else dict(c=c, N=N, primes=primes, ticks=[float(t) for t in t0], eig=[float(mp.log10(abs(x))) for x in e0], sens={})
for p in (ONLY or primes):
    d = arb(DELTA); v, e = vacuum([arb(q) * (1 + d) if q == p else arb(q) for q in primes]); ov = abs(sum(a * b for a, b in zip(v, v0))); f = xi_hat(v)
    S = [(mp.findroot(f, t, solver="secant") - t) / mp.mpf(DELTA) for t in t0]; out["sens"][str(p)] = [float(mp.log10(abs(x))) for x in S]
    print(f"   prime {p} (1 - overlap = {mp.nstr(1 - ov, 2)}): log10 |d tick / d delta| = " + ", ".join(f"{float(mp.log10(abs(x))):+.1f}" for x in S), flush=True)
    json.dump(out, open(fn, "w"))
