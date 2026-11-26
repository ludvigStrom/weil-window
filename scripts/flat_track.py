"""The lowest flat (edge-blind) ordinary state of the odd half, as a function of the window c (notes/ladder.md, step 4).
The form is computed in ball arithmetic and rounded to doubles; ordinary eigenvalues (above 0.2) and their vectors are then
accurate in double precision (the tiny ones are noise and are ignored).
For each ordinary state below 1.3: eigenvalue, share of weight within 0.02 of the window's ends, band of scaling frequencies
(10%, 50%, 90% points of its weight). Prediction to test: the lowest flat state has its band just above W = 2 pi c, width about
sqrt(8 pi^2 c), and level about 2/sqrt(c).
Usage: python3 flat_track.py N c1 c2 ...   (appends to flat_track.json)"""
import sys, math, json, os, fcntl, time
import numpy as np
from flint import arb
import hp_block as hb
N = int(sys.argv[1])
for c in [float(x) for x in sys.argv[2:]]:
    t0 = time.time(); F = hb.Form(c); L = F.L
    nodes = hb.quad_nodes(L, [arb(q).log() for q, p in F.pp if q < c], 2 * math.pi * N / float(L.mid()))
    Hv = [F.H(s) for s, w in nodes]; S_, C2_ = hb.transforms(nodes, Hv, L, N); Qo = hb.gram(S_, C2_, F.H(hb.ZERO), F.H(L), L, N)[1]
    Q = np.array([[float(Qo[i, j].mid()) for j in range(N)] for i in range(N)]); Q = (Q + Q.T) / 2; E, V = np.linalg.eigh(Q)
    Lf = float(L.mid()); om = 2 * math.pi * (np.arange(N) + 1) / Lf; x = np.linspace(0, Lf / 2, 3001); S = math.sqrt(2 / Lf) * np.sin(np.outer(x, om)); dep = Lf / 2 - x
    rows = []
    THR = float(os.environ.get("THR", "0.2"))
    for j in np.where((E > THR) & (E < 1.3))[0][:14]:
        f = S @ V[:, j]; w = f * f; tot = np.trapz(w, x); edge = float(np.trapz(np.where(dep < 0.02, w, 0), x) / tot)
        cw = np.cumsum(V[:, j] ** 2); band = [float(om[np.searchsorted(cw, q)]) for q in (0.1, 0.5, 0.9)]
        rows.append(dict(eps=float(E[j]), edge=edge, band=band))
    rec = dict(c=c, N=N, L=Lf, W=2 * math.pi * c, ntiny=int((E < THR).sum()), states=rows)
    with open("flat_track.lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX); R = json.load(open("flat_track.json")) if os.path.exists("flat_track.json") else []; R.append(rec); json.dump(R, open("flat_track.json", "w"))
    print(f"c = {c}: {time.time() - t0:.0f} s, {rec['ntiny']} tiny, lowest ordinary {rows[0]['eps']:.4f}", flush=True)
