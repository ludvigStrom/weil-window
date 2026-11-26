# The solved case: curves over finite fields (reading, 26 Nov evening)

Question (Ludvig, after "what would Connes do"): does our picture of the
window form exist where the Riemann hypothesis is a theorem, and if so,
what in it is the proved positivity? Read through a summarising fetch.

## What was read

**Suzuki, "A spectral theory of zeta functions of function fields",
arXiv:2610.41115 (2 Nov 2026).** The window form for a curve C over F_q
of genus g >= 1, already built.
- With h = log q and normalised Frobenius eigenvalues alpha_nu on the
  unit circle, p_n = sum_nu alpha_nu^-n; p_0 = 2g and for n >= 1
  p_n = (q^n + 1 - N_n) / q^(n/2), N_n = #C(F_{q^n}).
- The form on a window (-a, a) reduces to Toeplitz matrices
  H_r = (p_{l-k}), 0 <= k, l <= r, with r about 2a/h. The lowest
  eigenvalue is a step function of a, positive up to a finite a* and
  exactly zero beyond (there are only finitely many zeros).
- **Positivity is proved without using the location of the zeros**: with
  F the normalised Frobenius on the Jacobian, dagger the Rosati
  involution (F^dagger = F^-1) and tau the trace on the Tate module,
  tau(F^n) = p_n (Lefschetz) and tau(X X^dagger) >= 0 with equality only
  at X = 0 (Swinnerton-Dyer, by intersection theory on C x C). So
  sum a_k conj(a_l) p_{k-l} = tau(X X^dagger) for X = sum a_k F^k: the
  window form is a Gram matrix.
- Also: self-adjoint extensions of i d/dx on the window whose eigenvalues
  decrease to the zeros as a -> a*; translation by log q is unitarily
  equivalent to multiplication by Frobenius. Example: y^2 + y = x^3 + x
  + 1 over F_2: N_1 = 1, p_0 = 2, p_1 = sqrt 2; lowest eigenvalue 2h,
  (2 - sqrt 2) h, then 0.
- Not there: leaving out a point; genus 0.

**Burnol, arXiv:math/9911175.** A stationary process with spectral
density |Z(z/sqrt q)|^2; Lax-Phillips scattering; the scattering is
causal if and only if RH holds for the curve. **Causality is not proved
independently there**; it follows from Weil's theorem.

**Connes-Consani-Marcolli, arXiv:math/0703392.** A "very tentative"
dictionary between Weil's proof and the adele class space: Frobenius
correspondence <-> scaling correspondence; trivial correspondences <->
the range of a restriction map; Riemann-Roch <-> an index theorem;
positivity Tr(Z * Z') >= 2 d'(Z) >= 0 from Riemann-Roch on the curve.
Missing for number fields, by their account: the analogue of principal
divisors, a surface C x C, intersection theory.

## What this says about our picture

1. **The skeleton carries over.** A window of r + 1 sites has r + 1 modes;
   the curve has 2g zeros; H_r is the Toeplitz matrix of 2g point masses,
   so it has rank at most 2g and r + 1 - 2g exactly null states: "modes
   minus zeros", with exact zeros where we have tiny eigenvalues. Below
   2g sites the lowest eigenvalue is how well a function on the window
   can avoid the zeros (ladder.md).
2. **The pole and the points.** p_n is (pole terms q^(n/2) + q^(-n/2))
   minus (points N_n q^(-n/2)). For the projective line N_n = q^n + 1
   and the form is identically zero: the points cancel the pole exactly.
   Leaving out one point of degree 1 adds q^(-|n|/2) for n != 0 and
   nothing on the diagonal: a Toeplitz form with symbol (Poisson kernel
   minus 1), which is indefinite. So the analogue of "leave out a prime
   and the form fails" exists and is exactly computable. Not yet done.
3. **What the proved positivity is.** The form is a sum of squares
   because the clock (Frobenius) is unitary for a positive inner product
   on a finite-dimensional space (the Jacobian's H^1), and that inner
   product comes from geometry (intersection numbers on C x C). In the
   number field case the same sentence would need a space on which
   scaling is unitary with a positive inner product: the missing
   cohomology. Andrews-Connes-La Gatta's Tr(theta(f)(1 - P - Phat)
   theta(f)^*) has the Gram shape but 1 - P - Phat is not positive.
4. **A sobering point.** In the solved case the scattering picture
   (Burnol) is still only equivalent to RH; it does not prove it. The
   proof needed the surface. So a physical picture of our kind should not
   be expected to prove positivity even when it is true; what it can do
   is say what the missing positive object must look like.

## Next

Compute the toy: Toeplitz forms for an elliptic curve and a genus-2 curve
over a small field; the states (null, ordinary) against the zeros; leave
out a point and find when and how the form fails; compare with the
margin, the diver and the echoes of the number-field notes.

## The toy, computed

`scripts/track1_test_functions/ff_toy.py p "coefficients of f"` (output
`ff_toy.out`): curves y^2 = f(x) of genus 1 to 4 over F_2, F_3, F_5,
points counted by brute force over F_{p^n}, n <= g; zeta from Newton's
identities and the functional equation. `diver_floor.py`: the same
question for the number-field form.

**The curves behave as Suzuki says.** All zeros on the circle; the window
form on r + 1 sites is positive for r + 1 <= 2g and has an exact null
space beyond. Example, genus 4 over F_3 (y^2 = x^9 + x^4 + x + 2,
N_1..N_4 = 3, 15, 18, 75): lowest eigenvalue 8, 7.42, 6.00, 4.41, 3.79,
2.54, 1.36, 0.90 for 1 to 8 sites, then 0.

**Leaving out one point.** A closed point of degree d contributes
-d q^(-n/2) to p_n for d | n, n >= 1. Without it:

| curve | 2g | first window (sites) at which the form is negative, without one point of degree 1 / 2 / 3 | lowest eigenvalue at the largest window computed |
|---|---|---|---|
| genus 1 over F_2 | 2 | 4 / 5 / 4 | -0.821 / -0.921 / -1.255 |
| genus 2 over F_3 | 4 | 6 / 7 / 5 | -0.726 / -0.859 / -0.763 |
| genus 2 over F_5 | 4 | 4 / 5 / 4 | -0.610 / -0.428 / -0.408 |
| genus 3 over F_3 | 6 | 7 / 9 / 7 | -0.728 / -0.822 / -0.829 |
| genus 4 over F_3 | 8 | 9 / 10 / 9 | -0.729 / -0.622 / -0.856 |

- The form fails only once the window has null states, at 2g + 1 sites
  or one or two sites later. While the form is positive definite
  (at most 2g sites) leaving out a point shifts the lowest eigenvalue by
  a modest amount, up or down.
- **The fall is bounded, and the bound is a property of the removed term
  alone.** The removed term is a Toeplitz form with symbol
  d (P_r(d theta) - 1), P_r the Poisson kernel, r = q^(-d/2); its minimum
  is -2 d r/(1 + r): -0.732 for a degree-1 point over F_3 (computed:
  -0.726 to -0.729), -0.828 over F_2 (-0.821), -0.618 over F_5 (-0.610).
  Reason: beyond 2g sites the null space is large enough to contain a
  state concentrated where the removed symbol is most negative.

**The same in the number-field form.** Leaving out all powers of a prime
p (`diver_floor.py 200 p c...`, both halves):

| prime | c = p + 0.3 | larger windows | predicted long-run floor -2 ln p/(sqrt p + 1) |
|---|---|---|---|
| 2 | -0.12 | -0.49 (c = 5), -0.547 (10), -0.557 (20), -0.566 (50) | -0.574 |
| 3 | -0.14 | -0.63 (6), -0.65 (10), -0.73 (20), -0.76 (50) | -0.804 |
| 7 | -0.12 | -0.67 (10), -0.735 (14), -0.7355 (20, 30), -0.741 (50) | -1.068 |

- For 7 the level sits at -0.7355 = ln 7 / sqrt 7 from c = 14 until
  c = 49. While only the first power of p is in the window the removed
  term is lambda_p (T + T*) with T^2 = 0, whose lowest eigenvalue is
  exactly -lambda_p. When 49 enters it starts moving toward the
  long-run floor.
- So the "diver" of ladder.md and margin.md does not fall without limit.
  It falls from zero to minus the prime's own weight within about one
  unit of window, and stays there. The tiny states can be combined into
  a state on which the prime's mirror has eigenvalue exactly -1 (or +1,
  by parity): a state that is its own mirror image up to sign.

**What the toy says about our three objects.**
- *Margin.* In the toy there is no continuous margin, only "zero, one or
  two sites after null states appear". The continuous margin of the
  number-field form comes from its tiny eigenvalues being small but not
  zero.
- *Diver.* Bounded by the removed term's own symbol; the full form's
  near-null space is what lets the bound be reached.
- *Positivity.* In the toy the full form is a Gram matrix (Rosati), and
  each point's term is not positive by itself: only the complete sum
  over points, with the pole terms and the genus on the diagonal, is.
  Nothing in the toy suggests positivity can be assembled point by
  point, which matches induction.md.
