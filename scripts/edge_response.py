"""The edge response rho(w) of the odd half (notes/margin.md, part 5).
Window c, L = log c, shift s = L - w. On odd f the overlap S(f) = int f(y) f(y + s) dy only involves the two edge strips of
width w.  rho(w) = sup over odd f of -S(f) / Q(f) = top eigenvalue of Qo^-1 M,  M = -S.   A term of weight lambda at lag s,
if left out, makes the form fail exactly when 2 lambda rho = 1.
Everything in ball arithmetic (flint). Power iteration on X = Qo^-1 M.
Usage: python3 edge_response.py N c w1 [w2 ...]        (add CHECK=n to compare M with the withheld-term matrix of hp_block)"""
import sys, math, os, time, json
from flint import arb, arb_mat
import hp_block as hb
N = int(sys.argv[1]); c = float(sys.argv[2]); ws = [float(x) for x in sys.argv[3:]]
def odd_form(excl=()):
    F = hb.Form(c); L = F.L
    nodes = hb.quad_nodes(L, [arb(q).log() for q, p in F.pp if q < c], 2 * math.pi * N / float(L.mid()))
    Hv = [F.H(s, exclude=excl) for s, w in nodes]; S_, C2_ = hb.transforms(nodes, Hv, L, N)
    Qe, Qo = hb.gram(S_, C2_, F.H(hb.ZERO, exclude=excl), F.H(L, exclude=excl), L, N); return Qo, L
def Smat(L, s):
    n = N; l = L / 2; w = L - s; y0 = -l; y1 = -l + w; om = [2 * arb.pi() * (j + 1) / L for j in range(n)]
    rows = [[None] * n for _ in range(n)]
    for j in range(n):
        for k in range(n):
            a, b = om[j], om[k]
            f1 = w * (b * s).cos() if j == k else (((a - b) * y1 - b * s).sin() - ((a - b) * y0 - b * s).sin()) / (a - b)
            f2 = (((a + b) * y1 + b * s).sin() - ((a + b) * y0 + b * s).sin()) / (a + b)
            rows[j][k] = (f1 - f2) / L
    return arb_mat(n, n, [(rows[j][k] + rows[k][j]) / 2 for j in range(n) for k in range(n)])
t0 = time.time(); Qo, L = odd_form()
if os.environ.get("CHECK"):
    nn = int(os.environ["CHECK"]); Qw, _ = odd_form((nn,)); lam = math.log({4: 2, 8: 2, 9: 3}.get(nn, nn)) / math.sqrt(nn); S = Smat(L, arb(nn).log())
    D = Qw - Qo; err = max(abs(float((D[i, j] - 2 * lam * S[i, j]).mid())) for i in range(N) for j in range(N)); big = max(abs(float(D[i, j].mid())) for i in range(N) for j in range(N))
    print(f"check: largest entry of the withheld-term matrix {big:.3e}; largest difference from 2 lambda S: {err:.3e}")
out = []
for w in ws:
    M = -Smat(L, L - arb(repr(w))); X = Qo.solve(M)
    v = arb_mat(N, 1, [arb(1) / (i + 1) for i in range(N)]); rho = None; sig = float(os.environ.get("SHIFT", "0"))
    for it in range(int(os.environ.get("ITMAX", "400"))):
        u = X * v + v * arb(repr(sig)); nrm = max(abs(float(u[i, 0].mid())) for i in range(N)); i0 = max(range(N), key=lambda i: abs(float(u[i, 0].mid())))
        new = float((u[i0, 0] / v[i0, 0]).mid()) - sig if abs(float(v[i0, 0].mid())) > 0 else 0.0
        v = arb_mat(N, 1, [(u[i, 0] / arb(repr(nrm))).mid() for i in range(N)])   # midpoints only: radii would grow without bound over many iterations
        if rho is not None and abs(new - rho) < 1e-9 * abs(new): rho = new; break
        rho = new
    out.append((w, rho, it)); print(f"N = {N}, c = {c}, w = {w:.6f}: rho = {rho:.6f}  ({it} iterations, {time.time() - t0:.0f} s)", flush=True)
f = "edge_response.json"; R = json.load(open(f)) if os.path.exists(f) else []
R += [dict(N=N, c=c, w=w, rho=r) for w, r, _ in out]; json.dump(R, open(f, "w"), indent=0)
