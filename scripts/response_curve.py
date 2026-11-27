"""The response of the first tick to a mirror of unit strength at any distance (notes/questions.md, 'Extending the formula').
Same machinery as sens_fast.py. For a shift a = ln q (q need not be a prime power):
   T(q)  = d(first tick) for the perturbation  -A'_a   (moving a unit mirror),
   TA(q) = d(first tick) for the perturbation  -A_a    (strengthening it).
Then for any prime:  d tick / d delta = sum_k p^(-k/2) [ k ln p T(p^k) + (1 - (k/2) ln p) TA(p^k) ].
Usage: python3 response_curve.py N c q1,q2,...   (writes response_curve_<c>.json)"""
import sys, math, os, json, time
import mpmath
from mpmath import mp
from flint import arb
import hp_block as hb
mp.dps = int(os.environ.get("DPS", 170)); N = int(sys.argv[1]); c = float(sys.argv[2]); T0 = time.time()
F = hb.Form(c); L = F.L; tomp = lambda x: mp.mpf(x.mid().str(hb.DIG, radius=False)); Lm = tomp(L); n1 = N + 1
om = [2 * mp.pi * m / Lm for m in range(2 * N + 2)]; nrm = [1 / mp.sqrt(Lm)] + [mp.sqrt(2 / Lm)] * N
QS = [float(x) for x in sys.argv[3].split(",")]
nodes = hb.quad_nodes(L, [arb(q).log() for q, pp in F.pp], 2 * math.pi * N / float(L.mid()))
S_, C2_ = hb.transforms(nodes, [F.H(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, F.H(hb.ZERO), F.H(L), L, N)
Q = mp.matrix(n1, n1)
for i in range(n1):
    for j in range(n1): Q[i, j] = (tomp(Qe[i, j]) + tomp(Qe[j, i])) / 2
print(f"c = {c}, N = {N}: form built ({time.time() - T0:.0f} s)", flush=True); Qi = mp.inverse(Q); v = mp.matrix([mp.mpf(1) / (i + 1) for i in range(n1)])
for it in range(6): v = Qi * v; v = v / mp.norm(v)
lam0 = (v.T * Q * v)[0]; w2 = Qi * v; print(f"   vacuum eigenvalue {mp.nstr(lam0, 4)} (check: 1 / |Q^-1 v| = {mp.nstr(1 / mp.norm(w2), 4)}); inverse done ({time.time() - T0:.0f} s)", flush=True)
def hat(vec):
    def f(s):
        return sum(vec[m] * nrm[m] * (mp.sincpi((s - om[m]) * Lm / (2 * mp.pi)) + mp.sincpi((s + om[m]) * Lm / (2 * mp.pi))) for m in range(n1))
    return f
f0 = hat(v); GUESS = [mp.mpf(x) for x in os.environ.get("TICKS", "14.1347").split(",")]; ts = [mp.findroot(f0, g, solver="secant") for g in GUESS]; ds = [mp.diff(f0, t) for t in ts]
print("   ticks: " + ", ".join(mp.nstr(t, 12) for t in ts), flush=True)
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
if os.environ.get("DECOMP"):
    # the tilt: decompose the response over the near-null states above the vacuum.  TA_n(q) = sum_j c_j(q) m_jn,
    # c_j(q) = <v_j | A_a | v0>,  m_jn = v_j^(t_n) / (lambda_j v0^'(t_n))
    J = int(os.environ["DECOMP"]); vs = [v]; lams = [lam0]
    for j in range(1, J + 1):
        w = mp.matrix([mp.mpf(1) / (i + 2) ** 2 for i in range(n1)])
        for it in range(8):
            for u in vs: w = w - u * (u.T * w)[0]
            w = Qi * w
            for u in vs: w = w - u * (u.T * w)[0]
            w = w / mp.norm(w)
        vs.append(w); lams.append((w.T * Q * w)[0])
    print("   eigenvalues: " + ", ".join(mp.nstr(x, 3) for x in lams)); hats = [hat(u) for u in vs]
    m = [[hats[j](t) / (lams[j] * d) for t, d in zip(ts, ds)] for j in range(1, J + 1)]
    print("   m_jn relative to j = 1 (log10), per tick: " + " | ".join(f"j={j + 1}: " + " ".join(f"{float(mp.log10(abs(m[j][n] / m[0][n]))):+.1f}" for n in range(len(ts))) for j in range(1, J)))
    print("   q     share of each state j = 1.. in TA for tick 1            ... and for the last tick                 c_j(q)/c_1(q) (log10)")
    for q in QS:
        A, dA = Aops(mp.log(mp.mpf(repr(q)))); Av = A * v; cj = [(u.T * Av)[0] for u in vs[1:]]; row = []
        for n in (0, len(ts) - 1):
            terms = [cj[j] * m[j][n] for j in range(J)]; tot = sum(terms); big = max(range(J), key=lambda j: abs(terms[j])); row.append(f"sum 10^{float(mp.log10(abs(tot))):+.2f}, largest term j={big + 1} ({float(terms[big] / tot):+.2f} of the sum); " + " ".join(f"{float(x / tot):+.2f}" for x in terms[:6]))
        print(f"   {q:5.2f}  {row[0]}     {row[1]}     " + " ".join(f"{float(mp.sign(x / cj[0]) * 1):+.0f}e{float(mp.log10(abs(x / cj[0]))):+.1f}" for x in cj[1:]), flush=True)
    sys.exit()
def resp(vec):
    r = vec - v * (v.T * vec)[0]; dv = -(Qi * r); dv = dv - v * (v.T * dv)[0]; g = hat(dv); return [-g(t) / d for t, d in zip(ts, ds)]
if os.environ.get("OSC"):
    # test of the explanation: TA_n(a) should satisfy  TA'' + t_n^2 TA = c_n R(a),  R(a) = int v0(u) v0(u + a) du
    h = mp.mpf("1e-3"); print("   q      R(a) (log10)    (TA'' + t^2 TA) / R  for each tick")
    for q in QS:
        a0 = mp.log(mp.mpf(repr(q))); vals = []
        for da in (-h, 0, h):
            A, dA = Aops(a0 + da); vals.append(resp(-(A * v)))
            if da == 0: R = (v.T * A * v)[0] / 2
        rat = [((vals[0][n] - 2 * vals[1][n] + vals[2][n]) / h ** 2 + ts[n] ** 2 * vals[1][n]) / R for n in range(len(ts))]
        print(f"   {q:5.2f}  {float(mp.log10(abs(R))):+8.2f}     " + "  ".join(mp.nstr(x, 6) for x in rat), flush=True)
    sys.exit()
fn = f"response_multi_{c:g}.json"; out = json.load(open(fn)) if os.path.exists(fn) else {}; lg = lambda x: (float(mp.sign(x)), float(mp.log10(abs(x))))
print("   q      log10|TA| for each tick (sign)            ratio TA_n / TA_1 (log10)          T / TA for tick 1")
for q in QS:
    A, dA = Aops(mp.log(mp.mpf(repr(q)))); T = resp(-(dA * v)); TA = resp(-(A * v)); out[repr(q)] = dict(TA=[lg(x) for x in TA], T=[lg(x) for x in T])
    print(f"   {q:5.2f}  " + " ".join(f"{'+' if x > 0 else '-'}{float(mp.log10(abs(x))):+7.2f}" for x in TA) + "    " + " ".join(f"{float(mp.log10(abs(x / TA[0]))):+6.2f}" for x in TA[1:]) + f"     {mp.nstr(T[0] / TA[0], 4)}", flush=True); json.dump(out, open(fn, "w"))
