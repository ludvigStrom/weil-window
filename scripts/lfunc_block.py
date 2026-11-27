"""Weil's form on a window for a real primitive Dirichlet L-function (notes/questions.md, 'Another number system').
Same machinery as hp_block.py (kernel function H, transforms, Gram matrices). What changes, from the explicit formula for
L(s, chi) with conductor Nc and parity a (chi(-1) = (-1)^a):
  * each prime power q carries chi(q) (primes dividing Nc drop out);
  * there is no pole: the term -x U of zeta's H (measure -dt) is replaced, through the trivial zeros, by
        a = 0 :  -U                                  (psi(x, chi) contains -log x - (1/2) log(1 - x^-2))
        a = 1 :  -x U - U + (x + 1) log((x + 1) / 2)  (psi(x, chi) contains artanh(1 / x))
    i.e. zeta's H plus e^(-U/2) int_1^x (x/t - 1) m(t) dt with m = 1 + 1/t (a = 0) or 1 + 1/(t (t + 1)) (a = 1);
  * the conductor adds log(Nc) g(0), i.e. -(1/2) log(Nc) (x - 1) inside the bracket.
Check: the zeros of the transform of the lowest state should be the zeros of L(1/2 + i t, chi)."""
import math
from flint import arb
import hp_block as hb
CHARS = {"chi4": (4, 1, {1: 1, 3: -1}), "chi3": (3, 1, {1: 1, 2: -1}), "chi5": (5, 0, {1: 1, 2: -1, 3: -1, 4: 1}), "chi8": (8, 0, {1: 1, 3: -1, 5: -1, 7: 1}), "chi7": (7, 1, {1: 1, 2: 1, 3: -1, 4: 1, 5: -1, 6: -1})}
class FormChi(hb.Form):
    def __init__(self, lam2, name):
        super().__init__(lam2); self.name = name; self.Nc, self.a, tab = CHARS[name]; self.chi = lambda n: tab.get(n % self.Nc, 0)
        self.pp = [(q, p) for q, p in self.pp if self.chi(q) != 0]; self.logN = arb(self.Nc).log()
    def H(self, U, exclude=()):
        x = U.exp(); xf = float(x.mid()); S = hb.ZERO
        for q, p in self.pp:
            if q <= xf and q not in exclude: S += self.chi(q) * self.Lam[q] * (x / q - 1)
        pole = -U if self.a == 0 else -x * U - U + (x + 1) * ((x + 1) / 2).log()
        G = S + pole + (1 + hb.EG) * x - self.logN / 2 * (x - 1); y = 1 / x; omy = -((-U).expm1())
        T = hb.ONE - 2 * arb(2).log() / 2 if U.is_zero() else hb.ONE - (1 + y) / (2 * y) * y.log1p() + omy / (2 * y) * omy.log()
        return (-U / 2).exp() * (G - hb.LOG2PI + T)
