"""The window form for curves over finite fields, and what happens when a point is left out (notes/function_field.md).
Curve y^2 = f(x) over F_p, f of degree 2g+1 (one point at infinity). Points are counted by brute force over F_{p^n}, n <= g,
the zeta numerator P(T) = prod (1 - alpha_i T) follows from Newton's identities and the functional equation.
With normalised eigenvalues alpha/sqrt(q):  p_0 = 2g,  p_n = (q^n + 1 - N_n)/q^(n/2) = sum over zeros of e^(-i n theta).
Window of r+1 sites: Toeplitz matrix H_r = (p_{k-l}).  (Suzuki, arXiv:2610.41115; positive by the Rosati involution.)
Closed points: N_n = sum_{d | n} d a_d.  Leaving out one closed point of degree d adds d q^(-n/2) to p_n for d | n, n >= 1.
Usage: python3 ff_toy.py p "coefficients of f, lowest degree first"      e.g.  python3 ff_toy.py 3 "1 1 0 0 0 1" """
import sys, itertools, math
import numpy as np
p = int(sys.argv[1]); f = [int(x) % p for x in sys.argv[2].split()]; g = (len(f) - 2) // 2; q = p
def field(n):
    """GF(p^n) as tuples mod an irreducible polynomial found by search; returns elements, add, mul."""
    if n == 1: return [(a,) for a in range(p)], (lambda a, b: ((a[0] + b[0]) % p,)), (lambda a, b: ((a[0] * b[0]) % p,))
    def polymulmod(a, b, m):
        res = [0] * (2 * n - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b): res[i + j] = (res[i + j] + x * y) % p
        for i in range(2 * n - 2, n - 1, -1):
            c = res[i]
            if c:
                for j in range(n + 1): res[i - n + j] = (res[i - n + j] - c * m[j]) % p
        return tuple(res[:n])
    els = list(itertools.product(range(p), repeat=n))
    for low in itertools.product(range(p), repeat=n):
        m = list(low) + [1]; x = tuple([0, 1] + [0] * (n - 2)) if n > 1 else None
        # irreducible iff x^(p^n) = x and x^(p^k) != x for proper divisors k of n
        def frob_pow(k):
            y = x
            for _ in range(k):
                z = tuple([1] + [0] * (n - 1))
                for _ in range(p): z = polymulmod(z, y, m)
                y = z
            return y
        if m[0] != 0 and frob_pow(n) == x and all(frob_pow(k) != x for k in range(1, n) if n % k == 0):
            return els, (lambda a, b: tuple((s + t) % p for s, t in zip(a, b))), (lambda a, b, m=m: polymulmod(a, b, m))
    raise RuntimeError("no irreducible polynomial")
def count(n):
    els, add, mul = field(n); zero = tuple([0] * n); squares = {}
    for y in els: s = mul(y, y); squares[s] = squares.get(s, 0) + 1
    tot = 1                                                        # the point at infinity
    for x in els:
        v = zero; xp = tuple([1] + [0] * (n - 1))
        for c in f:
            if c: v = add(v, tuple((c * t) % p for t in xp))
            xp = mul(xp, x)
        tot += squares.get(v, 0)
    return tot
N = [None] + [count(n) for n in range(1, g + 1)]; S = [None] + [q ** n + 1 - N[n] for n in range(1, g + 1)]
e = [1.0]                                                           # elementary symmetric functions from power sums
for k in range(1, g + 1): e.append(sum((-1) ** (i - 1) * e[k - i] * S[i] for i in range(1, k + 1)) / k)
coef = [(-1) ** k * e[k] for k in range(g + 1)]; coef += [q ** (k - g) * coef[2 * g - k] for k in range(g + 1, 2 * g + 1)]   # P(T) = sum coef_k T^k
alpha = 1 / np.roots(coef[::-1]); theta = np.sort(np.angle(alpha)); 
print(f"y^2 = f(x) over F_{p}, f = {f} (lowest degree first), genus {g}"); print("  point counts N_1..N_g:", N[1:], "  |alpha|/sqrt(q):", np.round(np.abs(alpha) / math.sqrt(q), 6))
print("  angles of the zeros / pi:", np.round(theta / math.pi, 4)); gaps = np.diff(np.concatenate([theta, [theta[0] + 2 * math.pi]])); print("  gaps between neighbouring zeros / pi:", np.round(gaps / math.pi, 4))
def pn(n):
    n = abs(n); return 2.0 * g if n == 0 else float(np.real(np.sum((alpha / math.sqrt(q)) ** (-n))))
Nn = lambda n: round(q ** n + 1 - pn(n) * q ** (n / 2)); a = {}
for d in range(1, 2 * g + 6): a[d] = (Nn(d) - sum(k * a[k] for k in range(1, d) if d % k == 0)) // d
print("  closed points by degree a_1..a_6:", [a[d] for d in range(1, 7)])
def H(r, drop=None):
    M = np.array([[pn(k - l) for l in range(r + 1)] for k in range(r + 1)])
    if drop:
        for k in range(r + 1):
            for l in range(r + 1):
                n = abs(k - l)
                if n >= 1 and n % drop == 0: M[k, l] += drop * q ** (-n / 2)
    return M
R = 2 * g + 6; print(f"  lowest eigenvalue of the window form, r + 1 sites (2g = {2 * g} zeros):")
print("   sites | full form | " + " | ".join(f"minus one point of degree {d}" for d in (1, 2, 3) if a[d] > 0))
for r in range(0, R):
    row = [np.linalg.eigvalsh(H(r))[0]] + [np.linalg.eigvalsh(H(r, d))[0] for d in (1, 2, 3) if a[d] > 0]
    print(f"   {r + 1:5d} | " + " | ".join(f"{v:+10.5f}" for v in row))
