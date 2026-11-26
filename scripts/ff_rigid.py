"""Rigidity in the function-field toy (notes/questions.md). Window of r+1 sites, Toeplitz H = (p_{k-l}),
p_n = sum over zeros of e^(-i n theta). The 'primes' are the closed points: a_d of degree d, and changing a_d by x
changes p_n by -d x q^(-n/2) for every multiple n >= 1 of d. Parameters: a_1 .. a_r.
Which changes keep the form positive? Flat directions and the softness spectrum, as in rigid.py.
Usage: python3 ff_rigid.py q r theta1/pi theta2/pi ...   (one angle per conjugate pair)"""
import sys, math
import numpy as np
import soft
q = float(sys.argv[1]); r = int(sys.argv[2]); th = np.array([float(x) * math.pi for x in sys.argv[3:]]); g = len(th); ang = np.concatenate([th, -th])
pn = lambda n, a=ang: float(np.sum(np.cos(n * a))); T = lambda p: np.array([[p(abs(k - l)) for l in range(r + 1)] for k in range(r + 1)])
H = T(pn); E, U = np.linalg.eigh(H); V = U[:, E < 1e-9]; print(f"q = {q:g}, genus {g}, {r + 1} sites: {V.shape[1]} null states, lowest non-null eigenvalue {E[E > 1e-9][0]:.3f}")
M = [T(lambda n, d=d: -d * q ** (-n / 2) if n >= 1 and n % d == 0 else 0.0) for d in range(1, r + 1)]; P = [[V.T @ m @ V] for m in M]
print("   one degree alone (cost per unit change of a_d, up / down): " + ", ".join(f"d={d + 1}: {soft.lam(P, np.eye(r)[d]):+.3f}/{soft.lam(P, -np.eye(r)[d]):+.3f}" for d in range(r)))
s, F = soft.flat(P); print(f"   singular values of d -> sum d_i P_i: " + ", ".join(f"{x:.2e}" for x in s)); print(f"   exactly flat directions: {len(F)}   (2g - 1 = {2 * g - 1})")
# the known flat family: move the angles, or shift weight between zeros (total weight p_0 fixed)
tang = []
for j in range(g):
    dp = [0.0] + [-2 * n * math.sin(n * th[j]) for n in range(1, r + 1)]; tang.append(dp)                    # d p_n / d theta_j
    if j: tang.append([0.0] + [2 * math.cos(n * th[j]) - 2 * math.cos(n * th[0]) for n in range(1, r + 1)])   # weight from pair 0 to pair j
# convert a change of p_n into a change of a_d:  dp_n = -q^(-n/2) sum_{d|n} d da_d
def to_a(dp):
    da = {}
    for n in range(1, r + 1): da[n] = (-dp[n] * q ** (n / 2) - sum(d * da[d] for d in range(1, n) if n % d == 0)) / n
    return np.array([da[n] for n in range(1, r + 1)])
if len(F):
    Pr = F.T @ F; print("   'slide a zero' / 'shift weight between zeros', share lying in the flat space: " + ", ".join(f"{(lambda v: v @ Pr @ v / (v @ v))(to_a(t)):.6f}" for t in tang))
rng = np.random.default_rng(2); sp = soft.spectrum(P, rng, skip=list(F))
print("   softness of the remaining directions: " + ", ".join(f"{v:+.3f}" for v, d in sp))
for name, d in ([("flat direction 1", F[0])] if len(F) else []) + [("softest non-flat direction", sp[0][1])]:
    for size in (1e-3, 1e-2):
        da = size * d; lo = np.linalg.eigvalsh(H + sum(x * m for x, m in zip(da, M)))[0]; print(f"   true form along {name}, size {size:g}: lowest eigenvalue {lo:+.3e}")
