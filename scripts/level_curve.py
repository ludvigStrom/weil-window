"""The one curve: the vacuum's level lambda_0(c), even and odd halves, on a fine grid of windows (notes/questions.md).
Inverse iteration in mpmath on the ball-arithmetic form; appends to level_curve.json (resumable: skips windows already done).
Usage: python3 level_curve.py N c_start c_stop step     env DPS (default 70), HP_PREC, HP_DIGITS, CHI (chi4, chi3, ...: lfunc_block), LC_OUT, LC_TMAX"""
import sys, math, os, json, time
import mpmath
from mpmath import mp
from flint import arb
import hp_block as hb
CHI = os.environ.get("CHI")
if CHI: import lfunc_block as lf
mp.dps = int(os.environ.get("DPS", 70)); N = int(sys.argv[1]); c0, c1, st = (float(x) for x in sys.argv[2:5]); tomp = lambda x: mp.mpf(x.mid().str(hb.DIG, radius=False))
fn = os.environ.get("LC_OUT", "level_curve.json"); out = json.load(open(fn)) if os.path.exists(fn) else {}; T0 = time.time(); k = 0
def lowest(Q, n):
    M = mp.matrix(n, n)
    for i in range(n):
        for j in range(n): M[i, j] = (tomp(Q[i, j]) + tomp(Q[j, i])) / 2
    Mi = mp.inverse(M); v = mp.matrix([mp.mpf(1) / (i + 1) for i in range(n)])
    for it in range(5): v = Mi * v; v = v / mp.norm(v)
    return 1 / mp.norm(Mi * v)
while c0 + k * st <= c1 + 1e-9 and time.time() - T0 < float(os.environ.get("LC_TMAX", 150)):
    c = round(c0 + k * st, 6); k += 1; key = f"{c:.4f}|{N}" + (f"|{CHI}" if CHI else "")
    if key in out: continue
    F = lf.FormChi(c, CHI) if CHI else hb.Form(c); L = F.L; nodes = hb.quad_nodes(L, [arb(q).log() for q, pp in F.pp], 2 * math.pi * N / float(L.mid()))
    S_, C2_ = hb.transforms(nodes, [F.H(s) for s, w in nodes], L, N); Qe, Qo = hb.gram(S_, C2_, F.H(hb.ZERO), F.H(L), L, N)
    le, lo = lowest(Qe, N + 1), lowest(Qo, N); out[key] = [float(mp.log10(le)), float(mp.log10(lo))]; json.dump(out, open(fn, "w"))
    print(f"c = {c:.3f}: log10 lambda0 even {float(mp.log10(le)):.5f}, odd {float(mp.log10(lo)):.5f}   ({time.time() - T0:.0f} s)", flush=True)
