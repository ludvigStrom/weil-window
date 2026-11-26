"""Where the factor 4 to 5 of part 7 comes from (notes/margin.md, part 8).
For the K lowest odd states u_j (block inverse iteration in ball arithmetic):
  beta_j^2 = -(1/2) d eps_j/da from two neighbouring windows (Chekmachev's Hadamard formula),
  T_ij(w)  = <u_i, M(w) u_j>  exactly (M = minus the strip overlap),
  rho_K    = top eigenvalue of eps^-1/2 T eps^-1/2: the edge response carried by these K states, no thin-strip assumption.
Compared with the thin-strip formula beta^2 I_kappa(w), I_kappa(w) = int_0^w dt / sqrt(ln(1/(kappa t)) ln(1/(kappa (w - t)))),
for kappa = 1 (part 7) and kappa = 2 pi e^gamma (half-line calculation: symbol ln(|xi| / 2 pi), Mellin symbol 1/(2s) - gamma).
Usage: python3 edge_factor.py N c K w1 w2 ..."""
import sys, math, json, time
import numpy as np, mpmath
from mpmath import mp
from flint import arb, arb_mat
import hp_block as hb
mp.prec = hb.PREC
N = int(sys.argv[1]); c0 = float(sys.argv[2]); K = int(sys.argv[3]); ws = [float(x) for x in sys.argv[4:]]; t0 = time.time()
def form(c):
    F = hb.Form(c); L = F.L
    nodes = hb.quad_nodes(L, [arb(q).log() for q, p in F.pp if q < c], 2 * math.pi * N / float(L.mid()))
    Hv = [F.H(s) for s, w in nodes]; S_, C2_ = hb.transforms(nodes, Hv, L, N); return hb.gram(S_, C2_, F.H(hb.ZERO), F.H(L), L, N)[1], L
def tomp(x): return mp.mpf(x.mid().str(hb.DIG, radius=False))
def lowest(Qo):
    V = arb_mat(N, K, [arb(math.sin(1.7 * (i + 1) * (j + 1)) / (i + 1)) for i in range(N) for j in range(K)])
    for it in range(14):
        V = Qo.solve(V); cols = hb.gram_schmidt([[V[i, j] for i in range(N)] for j in range(K)], N)
        V = arb_mat(N, K, [cols[j][i] for i in range(N) for j in range(K)])
    H = V.transpose() * Qo * V; Hm = mpmath.matrix(K, K)
    for i in range(K):
        for j in range(K): Hm[i, j] = (tomp(H[i, j]) + tomp(H[j, i])) / 2
    E, W = mpmath.eigsy(Hm); order = sorted(range(K), key=lambda i: E[i])
    R = arb_mat(K, K, [arb(mpmath.nstr(W[i, order[j]], hb.DIG)) for i in range(K) for j in range(K)])
    return [E[i] for i in order], V * R
def Smat(L, s):
    n = N; l = L / 2; w = L - s; y0 = -l; y1 = -l + w; om = [2 * arb.pi() * (j + 1) / L for j in range(n)]; rows = [[None] * n for _ in range(n)]
    for j in range(n):
        for k in range(n):
            a, b = om[j], om[k]
            f1 = w * (b * s).cos() if j == k else (((a - b) * y1 - b * s).sin() - ((a - b) * y0 - b * s).sin()) / (a - b)
            f2 = (((a + b) * y1 + b * s).sin() - ((a + b) * y0 + b * s).sin()) / (a + b); rows[j][k] = (f1 - f2) / L
    return arb_mat(n, n, [(rows[j][k] + rows[k][j]) / 2 for j in range(n) for k in range(n)])
dc = 0.004 * c0; Q0, L = form(c0); E0, U = lowest(Q0); Em, _ = lowest(form(c0 - dc)[0]); Ep, _ = lowest(form(c0 + dc)[0])
da = (math.log(c0 + dc) - math.log(c0 - dc)) / 2; beta2 = [float(-(Ep[j] - Em[j]) / da / 2) for j in range(K)]; eps = [float(e) for e in E0]; Lf = float(L.mid())
print(f"c = {c0}, N = {N}, K = {K}  ({time.time() - t0:.0f} s)"); print("  j | eps_j | beta_j^2 | beta_j^2/eps_j")
for j in range(K): print(f"  {j + 1} | {eps[j]:.4e} | {beta2[j]:.4e} | {beta2[j] / eps[j]:.2f}")
KAP = 2 * math.pi * math.exp(0.5772156649)
json.dump(dict(c=c0, N=N, eps=eps, beta2=beta2), open(f"edge_rates_{c0}.json", "w"))
def I(w, kap):
    t = (np.arange(400000) + 0.5) / 400000 * w; return float((1 / np.sqrt(np.log(1 / (kap * t)) * np.log(1 / (kap * (w - t))))).sum() * w / 400000)
# profile of each state near the right edge: u(l - t) * sqrt(ln(1/(kappa t))) / beta should be 1
Uf = np.array([[float(U[i, j].mid()) for j in range(K)] for i in range(N)]); om = 2 * math.pi * (np.arange(N) + 1) / Lf
print("  profile check: |u_j(edge - t)| / beta_j times sqrt(ln(1/t)) and times sqrt(ln(1/(kappa t))), kappa = 2 pi e^gamma = %.2f" % KAP)
for t in (0.003, 0.006, 0.012, 0.024, 0.048):
    if KAP * t >= 1: continue
    u = math.sqrt(2 / Lf) * np.sin(om * (Lf / 2 - t)) @ Uf
    print(f"     t = {t:.3f}: " + "  ".join(f"{abs(u[j]) / math.sqrt(beta2[j]) * math.sqrt(math.log(1 / t)):.2f}/{abs(u[j]) / math.sqrt(beta2[j]) * math.sqrt(math.log(1 / (KAP * t))):.2f}" for j in range(K)))
E = [r for r in json.load(open("edge_response.json")) if abs(r["c"] - c0) < 1e-6 and r["N"] == N]
print("  w | rho (full, this N) | rho_K (these K states) | trace sum T_jj/eps_j | formula kappa=1 | formula kappa=2 pi e^gamma | per state T_jj / (beta_j^2 I_kappa)")
for w in ws:
    M = -Smat(L, L - arb(repr(w))); T = U.transpose() * M * U; Tm = mpmath.matrix(K, K)
    for i in range(K):
        for j in range(K): Tm[i, j] = (tomp(T[i, j]) + tomp(T[j, i])) / 2 / mp.sqrt(E0[i] * E0[j])
    rK = float(max(mpmath.eigsy(Tm, eigvals_only=True))); tr = float(sum(Tm[j, j] for j in range(K)))
    full = next((r["rho"] for r in E if abs(r["w"] - w) < 1e-12), float("nan")); S = sum(beta2[j] / eps[j] for j in range(K))
    f1 = I(w, 1) * S; f2 = I(w, KAP) * S if KAP * w < 1 else float("nan")
    per = "  ".join(f"{float(tomp(T[j, j])) / (beta2[j] * I(w, KAP)):.2f}" for j in range(K)) if KAP * w < 1 else ""
    print(f"  {w:.5f} | {full:.3f} | {rK:.3f} | {tr:.3f} | {f1:.3f} | {f2:.3f} | {per}", flush=True)
print(f"  ({time.time() - t0:.0f} s)")
