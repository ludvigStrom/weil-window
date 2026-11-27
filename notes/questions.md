# Questions the physical picture suggests (started 26 Nov evening)

List agreed with Ludvig: (1) pressure of the gas, (2) entanglement
entropy, (3) predict a level from the zeros alone, (4) close pairs of
zeros, (5) self-mirrored states, (6) two mirrors at once, (7) heat
capacity, (8) what radiates, (9) force on a prime's position. Order of
attack: 3, then 2 and 6.

## Q3: levels from the zeros alone. Answer: yes, exactly, and cheaply

Weil's form is the sum over zeros. For odd real f on the window,
Q(f) = 2 sum over gamma > 0 of |int f(x) sin(gamma x) dx|^2. In the sine
basis Q = 2 A A^T with A[k, gamma] in closed form. So the ordinary
spectrum needs only a list of zeros and a few lines of numpy.

`level_from_zeros.py c N` (zeros below 2300 from the Riemann-Siegel
formula in doubles):

| c = 100.2: true form (ball arithmetic, N = 900) | from the zeros alone | state's median frequency |
|---|---|---|
| 0.1947 | 0.1944 | 651.9 |
| 0.2471 | 0.2472 | 542.8 |
| 0.2662 | 0.2663 | 725.5 |
| 0.2824 | 0.2808 | 631.4 |
| 0.3029 | 0.3028 | 897.4 |
| 0.3076 | 0.3074 | 563.2 |
| 0.3308 | 0.3306 | 788.3 |

Same at c = 50.2 (0.2694 / 0.2691, 0.3364 / 0.3360, ...). This is also an
independent check of `hp_block.py`: two unrelated computations of the
same form agree to three or four digits. (The tiny eigenvalues cannot be
got this way: they need the zeros exactly.)

**Local model.** Keep only the modes within 25 of a gap's centre and the
zeros within 120 of it: 0.2005 for the gap 650.67-653.65 (true 0.1947),
0.315 for 896.63-899.22 (0.303), 0.357 for 787.47-790.06 (0.331). So a
level is fixed to 3 to 8% by the zeros within about a hundred of its gap.

**The bottom of the ordinary spectrum for c = 20 to 160**
(`mu_from_zeros.py 2.0 c...`, zeros below 6000, one minute for 36
windows; data `mu_from_zeros.json`):

| windows c | gap in use | gap in mean spacings | level, from - to | level x sqrt(c)/2 |
|---|---|---|---|---|
| 20 | 107.2-111.0 | 1.75 | 0.44 | 0.99 |
| 24-28 | 169.9-173.4 | 1.84 | 0.39 - 0.37 | 0.96 - 0.99 |
| 32-36 | 279.2-282.5 | 1.96 | 0.40 - 0.39 | 1.14 - 1.17 |
| 40-68 | 415.5-418.4 | 1.96 | 0.33 - 0.23 | 1.05 - 0.91 - 0.95 |
| 72-108 | 650.7-653.6 | 2.20 | 0.24 - 0.194 - 0.20 | 1.03 - 0.97 - 1.05 |
| 112-132 | 724.6-727.1 | 1.88 | 0.22 - 0.24 | 1.15 - 1.40 |
| 136-148 | 1267.6-1270.1 | 2.15 | 0.24 - 0.21 | 1.41 - 1.28 |
| 152-160 | 1055-1070 | 1.4 to 2.1 | 0.20 - 0.17 | 1.22 - 1.05 |

- The level follows one gap for a stretch, improving as the window
  lengthens, then hands over to another. A gap goes on serving for a
  while after 2 pi c has passed it.
- 2/lambda holds to 10% from c = 20 to 110 and is 15 to 40% too low from
  112 to 150. Andrews-Connes-La Gatta's 0.19 to 0.20 at 150 is on our
  curve (0.21 at 148, 0.198 at 152).
- The gaps in use are 1.75 to 2.2 mean spacings wide.

## Q6: two mirrors at once. Answer: the floors add

`diver_floor.py N p1,p2 c...`: lowest eigenvalue with all powers of two
primes left out.

| primes out | sum of first-power floors (ln p / sqrt p) | sum of long-run floors | computed |
|---|---|---|---|
| 3 and 5 | -1.354 | -1.799 | -1.355 (even half, c = 24), -1.254 (odd) |
| 5 and 7 | -1.455 | -2.062 | -1.385 (odd, c = 24), -1.205 (even) |
| 2 and 3 | -1.124 | -1.378 | -1.252, -1.267 (c = 50, where 4 to 32 and 9, 27 are in) |

The lowest level heads for the sum of the two single floors: a state can
be its own image in both mirrors at once. No resonance beyond the sum
(none is possible: the removed terms' norms add). Consistent with the
mirrors being independent.

## Q2: entanglement entropy. Not well posed as asked

The entropy of the gas uses occupation numbers between 0 and 1, which
are the prolate concentration eigenvalues, not the Weil eigenvalues (the
ordinary Weil levels exceed 1 and grow like a logarithm). Computing it
from the prolate values would only reproduce Eisler-Peschel's known
result for a Fermi gas. Dropped unless a sensible map from Weil levels to
occupations turns up.

## Q3 continued: large windows (`mu_big.py`, zeros below 21000)

Same computation, full sine basis up to 2 x 2 pi c (checked unchanged with
4 x at c = 250), and at c = 1000 a band of modes 0.75 to 1.8 x 2 pi c
(a band reproduces the full answer at c = 500: 0.0556 against 0.0555).
"Lowest gap state" = lowest eigenvalue above 0.02 whose state sits above
2 pi c.

| c | lowest gap state | its gap (mean spacings) | at | 2/lambda | level x c |
|---|---|---|---|---|---|
| 100.2 | 0.194 | 650.67-653.65 (2.20) | 1.04 W | 0.200 | 19 |
| 160.2 | 0.166 | near 1060 | | 0.158 | 27 |
| 200.2 | 0.124 | 1331.83-1333.67 (1.57) | 1.06 W | 0.141 | 25 |
| 250.2 | 0.107 | 2817.95-2820.49 (2.47) | 1.79 W | 0.126 | 27 |
| 300.2 | 0.079 | 1975.08-1977.17 (1.92) | 1.05 W | 0.115 | 24 |
| 400.2 | 0.069 | 2817.95-2820.49 (2.47) | 1.12 W | 0.100 | 27 |
| 500.2 | 0.056 | 3824.54-3826.95 (2.46) | 1.22 W | 0.089 | 28 |
| 1000.2 | 0.030 | 7005.10-7006.74 (1.83) | 1.11 W | 0.063 | 30 |

- **The level keeps falling, and faster than 2/lambda.** From c = 160 to
  1000 it is closer to 27/c than to 2/sqrt(c). 2/lambda is 15 to 40% too
  low at 112-150, right near 160-200, and twice too high at 1000.
  Eight points only; "27/c" is a description, not a law.
- **There is no sharp line between ordinary and tiny at large c.** A
  state in a given gap falls steadily as the window grows (gap 1975-1977:
  0.110 at c = 250, 0.079 at 300, 0.058 at 400) and ends up among the
  tiny ones once 2 pi c has passed it. Between 0.02 and 0.1 there are a
  handful of states at every large window, some localised in a gap, some
  broad ("fallers"). So "the bottom of the ordinary spectrum" needs a
  definition before it can have an asymptotic law. How Andrews-Connes-
  La Gatta condition theirs is the thing to read next.
- Very wide gaps far above 2 pi c compete: at c = 250 the 2.47-spacing gap
  at 1.79 W wins. Higher up, zeros are denser relative to the window's
  resolution (ratio ln c / ln(gamma / 2 pi)), so a gap has to be wider to
  serve.

## Q9: the force on a prime's position. Answer: a sharp tent, peaked at the true prime

`prime_move.py N c p delta...`: the prime p is replaced by a "prime" at
x = p (1 + delta), weight ln x, powers x^k; nothing else changes. Lowest
eigenvalue of each half (double precision, so the true value at delta = 0,
far below 1e-15, reads as zero).

Prime 7, window c = 14.2 (N = 100; N = 50 and 200 give the same slopes):

| delta | position | lowest eigenvalue (even, odd) | per unit delta |
|---|---|---|---|
| -0.1 | 6.3 | -0.815, -0.760 | 8.2, 7.6 |
| -0.01 | 6.93 | -0.137, -0.132 | 13.7, 13.2 |
| -0.001 | 6.993 | -0.0134, -0.0135 | 13.4, 13.5 |
| 0 | 7 | 0 | |
| +0.00001 | 7.00007 | -0.000133, -0.000136 | 13.3, 13.6 |
| +0.001 | 7.007 | -0.0145, -0.0137 | 14.5, 13.7 |
| +0.01 | 7.07 | -0.141, -0.133 | 14.1, 13.3 |
| +0.03 | 7.21 | -0.430, -0.428 | 14.3, 14.3 |
| +0.1 | 7.7 | -0.642, -0.613 | 6.4, 6.1 |

- **The lowest eigenvalue is -k |delta|: negative on both sides, linear
  from the first step, a kink at the true position.** Positivity holds at
  one point only. It is linear, not quadratic, because the form has many
  near-zero states and the displacement acts on that whole near-null space
  with both signs: whichever way the prime moves, some state gains.
- The slope does not depend on the basis size, so it is a property of the
  window, carried by the tiny states (not by high frequencies).
- Saturates near the single-prime floor (-0.74 for 7) once the prime is
  several percent away: by then it is "7 missing and a stranger present".

Stiffness k (average of the two halves, right / left):

| window c | prime 2 | 3 | 7 | 11 | 13 | 23 |
|---|---|---|---|---|---|---|
| 7.2 | | 7.9 / 5.8 | 0.01 / 0.003 | | | |
| 14.2 | 47 / 37 | 25 / 27 | 14 / 13 | 13 / 10 | 7.7 / 5.2 | |
| 28.2 | | 91 / 90 | 53 / 39 | | | 24 / 18 |
| 56.2 | | | 123 / 115 | | | |

- Stiffer for small primes, and growing with the window faster than
  proportionally (prime 7: 14, 46, 119 at c = 14, 28, 56, roughly
  c^1.5; three points). No formula yet.
- A prime that has only just entered the window is hardly pinned at all
  (7 at c = 7.2: 0.01), the same delay as the margin in `string.md`.
- Physical reading: each prime sits at the tip of a tent, held from both
  sides with a force that does not vanish at the equilibrium (like a bead
  in a V-shaped groove, not a ball in a bowl), and a larger window holds
  it harder.
- Not checked against the literature. That positivity fails for wrong
  primes is expected by everyone; the tent shape and its stiffness may or
  may not have been computed.

## How Andrews-Connes-La Gatta define mu (read 26 Nov, via a summarising fetch of arXiv:2611.13129, section 5.3, eq. (40), Table 7)

mu^odd = minimum of the form over the odd Galerkin space, orthogonal to
lambda^2 - 2 prolate vectors and with the pole functional zero: lambda^2 - 1
conditions. So "ordinary" is defined **by count**, not by a threshold.
They note (Courant-Fischer) mu^odd <= the lambda^2-th odd eigenvalue of the
raw form; the (lambda^2 - 1)-th, "transitional", one is removed by the
pole condition.

Our raw spectra have exactly that shape: c - 2 eigenvalues below 0.02,
then one transitional (0.04 to 0.06), then the rest. So their mu^odd is
bounded above by our c-th eigenvalue:

| c | our c-th odd eigenvalue (from zeros) | their mu^odd (Table 7) | 2/lambda | lambda x value |
|---|---|---|---|---|
| 24 | 0.390 | 0.383 | 0.408 | |
| 28 | 0.374 | 0.361 | 0.378 | |
| 32 | 0.400 | 0.395 | 0.354 | |
| 36 | 0.387 | 0.377 | 0.333 | |
| 40 | 0.331 | 0.332 | 0.316 | |
| 50 | 0.269 | 0.270 | 0.283 | 1.91 |
| 100 | 0.1944 | 0.194 | 0.200 | 1.94 |
| 150 | 0.21 / 0.198 (148 / 152) | 0.194 | 0.163 | 2.4 |
| 200 | 0.124 | not computed | 0.141 | <= 1.75 |
| 250 | 0.0965 | | 0.126 | <= 1.53 |
| 300 | 0.079 | | 0.115 | <= 1.37 |
| 400 | 0.058 | | 0.100 | <= 1.15 |
| 500 | 0.057 | | 0.089 | <= 1.28 |

(Ours are at c + 0.2.) Where both exist they agree to a few percent, so
the bound is nearly tight there. Their stated range is 1.6 <= lambda mu
<= 2.5 for lambda^2 from 5 to 150. **Prediction: beyond 250 lambda mu
drops below 1.6 and keeps falling**; this follows from their own
inequality plus our eigenvalues, which need only the zeros. The earlier
worry (no sharp line between ordinary and tiny) is answered by their
count-based definition; note the c-th state at c >= 250 is often a state
in a gap *below* 2 pi c (250: gap 1330.4-1331.8 at 0.85 W), not above.

## Rigidity: moving all primes at once (`rigid.py`)

Every prime p <= c moved to p (1 + d_p). To first order the lowest
eigenvalue is lambda_min of sum d_p P_p on the near-null space V of the
true form (P_p = derivative of the form, restricted to V). We maximise
that over unit vectors d (concave, so local search suffices): if the
maximum is negative, no direction keeps the form positive.

c = 14.2, primes 2, 3, 5, 7, 11, 13, V = eigenvalues below 0.01 (13 even,
12 odd states):

- Each P_p alone has both signs (single-prime tents: 47/38, 27/28, 20/17,
  15/14, 14/11, 10/7 per unit displacement).
- A random direction: down to -53 per unit.
- **Best direction: -0.765 per unit. Still negative, so the true
  configuration is rigid, but 15 to 60 times softer than any single
  move.** Checked on the true form: -7.6e-4 at size 0.001, -7.6e-3 at
  0.01, exactly the first-order prediction, both halves.
- **The soft direction is d_p proportional to ln p**: found
  (0.158, 0.253, 0.373, 0.455, 0.529, 0.534), ln p normalised is
  (0.153, 0.243, 0.356, 0.430, 0.531, 0.567). That is p -> p^(1+eps):
  every log-position stretched by the same factor. It is the dilation of
  the log axis, the one near-symmetry of the form; it fails only because
  the archimedean part and the weights 1/sqrt(p) do not stretch with it
  (equivalently: the zeros would move to real part 1/(2(1+eps))).

With V = eigenvalues below 1e-7 only, the search instead returns "move
the newest prime alone" (-0.12), and the first-order number is wrong by
a factor 60 against the true form: the transitional states (1e-5 to
1e-2) must be in V. c = 7.2 (N = 80): same picture, the newest prime 7
is not yet held.

### Rigidity, continued: the softness spectrum, and c = 28.2

`rigid.py` now caches the derivatives (resumable) and, after the best
direction, finds the best direction orthogonal to it, and so on.

c = 14.2 (per unit displacement):

| direction | softness | what it is |
|---|---|---|
| 1 | -0.77 | all primes together, proportional to ln p (overlap 0.9989) |
| 2 | -7.9 | 5 up, 13 down |
| 3 | -11.6 | 7 down, 11 up |
| 4 | -14.7 | 5 down, 11 and 7 up |
| 5 | -17.8 | 3 alone |
| 6 | -39.1 | 2 alone |

**One soft direction, then a tenfold jump.** The others are neighbouring
primes moving against each other, hardest for the smallest primes.

c = 28.2, nine primes, N = 130, V = 27 + 26 states: zoom direction found
again, overlap with ln p 0.9998, softness -1.15 (true form at size 0.001:
-1.03e-3 and -1.14e-3; prediction -1.15e-3); single primes 20 to 122,
random directions down to -144. So from c = 14 to 28 the single-prime
stiffness triples and the zoom softness grows by half: the separation
widens with the window. (The orthogonal passes at c = 28 did not converge
in the time available; Nelder-Mead in 9 dimensions is the wrong tool,
use the subgradient or an SDP.)

The state that falls under the zoom is built mostly from the least-small
near-null states (the transitional ones, eigenvalues 1e-6 to 1e-14).

Reasoning, not computed: the form is 2AA^T >= 0 for *any* real positions
of the zeros, so sliding zeros along the line never breaks positivity.
A prime move breaks it because discrete primes cannot reproduce such a
slide. p -> p^(1+eps) is the one move of atoms that does slide all the
zeros (gamma -> gamma/(1+eps)); it also shifts their real part to
1/(2(1+eps)), and that shift, of size eps/2, is the small cost.

### Where the cost of the zoom comes from (`zoom_decomp.py`)

The guess "the cost is the zeros leaving the line by eps/2" was tested and
is **wrong in emphasis**. (It also fails on paper: for a real test
function a zero at real part 1/2 + eta contributes |F|^2 + O(eta^2); the
first-order term cancels between gamma and -gamma.)

Exact bookkeeping: stretching the test functions turns window c into
c' = c^(1/(1+eps)), where the true form is positive. So
Q_zoomed(c) = Q_true(c') + eps D, and D = D_smooth (the same difference
with all primes removed: archimedean and pole terms, which do not
stretch) + D_prime (prime weights: (1+eps) n^(-eps/2) instead of 1).

c = 14.2, on the falling state, per unit eps:

| | true form on c' | smooth terms | prime weights | total |
|---|---|---|---|---|
| zoom out, odd half (the one that falls first) | +0.01 | -3.62 | -0.06 | -3.67 |
| zoom out, even half | +0.00 | -1.00 | -0.84 | -1.84 |
| zoom in, even half (falls first) | +0.06 | -6.40 | +0.82 | -5.52 |
| zoom in, odd half (CORRECTED 26 Nov: this half does not fall; lowest eigenvalue is +0.57 per unit of shrink) | | | | +0.57 |

(Per unit eps; divide by |ln p vector| = 4.52 for "per unit displacement":
0.81 out, 1.22 in, against 0.77 found by the search.)

**The scale is held by the smooth background (the archimedean term and
the pole term), not by the primes.** The primes among themselves hold the
pattern; the place at infinity fixes the unit of length. The falling
state is a near-null state of the true form on the shrunk window, as it
should be.

## Rigidity in the curve setting, and the comparison (`ff_rigid.py`, `rigid_weights.py`, `soft.py`)

Toy: genus 2 curve y^2 = x^5 + x + 1 over F_3 (zeros at angles +-0.3041 pi,
+-0.6959 pi), window of r + 1 sites, Toeplitz form. Parameters: the numbers
a_d of closed points of degree d = 1..r (a change of a_d changes p_n for
all multiples n of d: a "prime" with its powers). Positions cannot move
(degrees are integers), so this is the analogue of changing weights.

| | 8 sites (4 null states) | 12 sites (8 null states) |
|---|---|---|
| exactly flat directions | 3 = 2g - 1 | 3 = 2g - 1 |
| what they are | slide a zero along the circle; shift weight between zeros (each lies 100% in the flat space) | same |
| directions that make the form strictly positive | yes (best +0.33 per unit) | yes (+0.05) |
| the rest | -0.13 to -0.62 | -0.05 to -1.16 |
| one degree alone | tent, 0.06 to 2.0 (except degree 4 downwards, which helps: +0.27) | tent, 0.02 to 2.3 |

So for a curve the true data sit on the boundary of a convex set with
interior: sliding zeros is free (other "curves"), smearing the zeros
(all p_n shrunk towards the trivial background p_0 = 2g) is strictly
allowed, and single points are pinned by tents just as primes are.

Number field, c = 14.2, same analysis (subgradient search; it reproduces
the earlier position spectrum):

| parameters | flat directions | softest | spectrum |
|---|---|---|---|
| weights only (6) | 0 | -0.15 (the newest prime's weight) | -0.15, -0.44, -0.59, -0.64, -0.72, -0.85 |
| positions only (6) | 0 | -0.77 (zoom) | -0.77, -7.9, -11.6, -14.7, -17.8, -39.1 |
| both (12) | 0 | -0.027 (zoom in, with weights lowered) | -0.027, -0.24, ..., -0.94, -7.9, ..., -39.1 |

- A single weight's tent has slope close to ln p / sqrt p (7: 0.74 / 0.73;
  5: 0.72 / 0.72): changing a weight is removing a fraction of the prime.
- All weights down together: -0.82. The background alone is not
  positive, so there is no "smearing" direction. In the curve case the
  background is 2g times the identity, which is.
- **No flat direction and no allowed direction at all**, to first order.
  (Strictly the form is positive definite, so there is an allowed ball,
  of radius about smallest eigenvalue / slope: 1e-50 here. That is the
  margin of `string.md`.)

**The difference between the two settings, as measured:** a curve has
neighbours (its zeros can slide, 2g - 1 free directions) and a positive
background to retreat to. The primes, within these parameters, have
neither: the nearest thing to a neighbour is the zoom, and the background
itself is negative without them.

Correction (26 Nov, late): the first run of `ff_rigid.py` used angles
0.2677 pi and 0.6177 pi, not the curve's 0.3041 pi and 0.6959 pi. The
table above is the rerun with the right angles; the structure (3 flat
directions, an allowed open set, tents) is the same for both.

Caveat: "these parameters" = positions and weights of the primes already
in the window. New primes, or a changed background, are not in the search.

## Literature pass before "other number systems" (26 Nov, search snippets only, nothing read in full)

- **Groskin, arXiv:2607.02828** (July 2026): every value of the truncated
  form is an exact, closed-form sum over the zeros (a "finite
  Guinand-Weil dictionary"). **So our Q = 2AA^T "from the zeros alone" is
  not new**; cite it in note v2. Whether he used it for mu or gaps: not
  seen.
- **Groskin, arXiv:2605.20224 and Zenodo 19546515**: independent
  implementation of the Connes-van Suijlekom matrix; a partial extension
  to L(s, chi_3) "shows convergence at c = 13, 17, 37 but an oscillating
  positivity breakdown at c = 23, 29". Someone has tried another
  L-function and hit trouble. Read before building ours.
- **Bombieri, Remarks on Weil's quadratic functional I (2000)**: if
  finitely many zeros are off the line, the number of negative eigenvalues
  is exactly half their number (truncation large enough). A test for our
  moved primes: count the negative eigenvalues.
- **Burnol, math/9809119, math/0101068**: explicit formula for any number
  field and Hecke character with all places on one footing; Weil
  positivity checked under a support condition.
- **Beurling primes**: Diamond-Montgomery-Vorhauer built a system where RH
  fails; Broucke-Vindas (arXiv:2004.11501) one where it holds. So other
  "prime systems" with zeros on a line exist, but with a different
  integer-counting side (background).
- Zenodo 20710075 (independent): prime contribution as a Loewner matrix,
  "the primes cooperate". Suzuki arXiv:2606.09096 (screw function) unifies
  Yoshida, Bombieri, Connes-Consani, CCM.
- Not found: anyone moving primes continuously and measuring the tent, the
  zoom direction, or a path between two number systems.

## Bombieri's test: how many levels go negative (`neg_count.py`, c = 14.2)

| move | negative eigenvalues, even / odd (of 13 / 12 near-null states) | sizes |
|---|---|---|
| prime 7 by +0.1% to +3% | 6 / 7 | spread, from -14 delta down to -1e-4 delta |
| prime 7 by +10%, or -1% | 6 / 6 | |
| prime 2 by +1% | 5 / 5 | |
| zoom out (p -> p^1.001, p^1.01) | 11 / 12 | all nearly equal: -1.4 to -1.8 eps (one odd: -3.7 eps) |
| zoom in (p -> p^0.99) | 1 / 0 | one even state, -5.7 eps |

- Same counts at N = 200 and for every size of move: the count is a
  property of the direction.
- **One prime moved: half of the near-null states fall** (about c/2 per
  half). In Bombieri's terms, as many zeros leave the line as the window
  can see.
- **Zoom out: all near-null states sink together, at one rate.** The
  smooth-term mismatch acts on the near-null space as (almost) a multiple
  of the identity.
- **Zoom in: everything lifts, except one even state.** The odd half
  stays positive. (This corrects the zoom table above, where a sign was
  misread.) One even direction is the natural place for the pole: the
  zoomed zeta has its pole at 1/(1+eps), the form has it at 1. If that
  one direction is the pole direction, then zooming in is allowed on the
  pole-free subspace. Not checked.

## Groskin's papers (summarising fetch; first 100k characters of 2605.20224, abstract page of 2607.02828)

- 2607.02828: Theorem 1 is the exact zero-sum representation of every
  value of the truncated form (even vectors). Theorem 2 certifies
  positivity against the archimedean cutoff. Nothing on mu, 2/lambda or
  gaps between zeros in what was read.
- 2605.20224: the chi_3 form works (16 to 28 digits at c = 13, 17, 37).
  The single negative even eigenvalue at c = 23, 29 is, by the author's
  own revision, an artifact of cutting the archimedean integral, gone at
  larger cutoff. Bulk eigenvalues 2..101 reported Poisson-like. How the
  form was modified for chi_3 was not in the part read.

## Zoom in, with the pole direction removed (`zoom_pole.py`)

Pole vectors as in ACL: 2 cosh(t/2) (even), 2 sinh(t/2) (odd).

| move | negative, even / odd | falling state's overlap with the pole vector | negative without the pole direction | lowest eigenvalue without it |
|---|---|---|---|---|
| c = 14.2, zoom in 1% | 1 / 0 | 0.977 (even) | 0 / 0 | +3.8e-3, +5.8e-3 |
| c = 14.2, zoom in 0.1% | 1 / 0 | 0.972 | 0 / 0 | +4.0e-4, +5.8e-4 |
| c = 28.2, zoom in 1% | 1 / 0 | 0.967 | 0 / 0 | +2.3e-3, +2.7e-3 |
| c = 14.2, zoom out 1% | 11 / 12 | 0.006 / 0.965 (odd: the pole falls too) | 11 / 11 | -0.018, -0.019 |
| c = 28.2, zoom out 1% | 25 / 25 | 0.003 / 0.983 | 25 / 24 | -0.021, -0.021 |
| c = 14.2, prime 7 by 1% | 6 / 7 | 0.04 / 0.07 | 6 / 6 | -0.14, -0.13 |

**Confirmed: the one state that falls when the primes are zoomed in is
the pole direction (97%). On the pole-free subspace, p -> p^(1-eps) keeps
the form positive, and makes it strictly positive: every near-null level
lifts to at least about 0.25 to 0.4 eps.** Zooming out does the opposite:
all near-null levels sink together at about 1.8 eps (c = 14) to 2.1 eps
(c = 28), and the pole direction switches halves (odd).

Reasoning, not computed: the archimedean term is logarithmic in frequency
(symbol about log |omega|), so stretching the frequency axis by 1 + eps
changes it by a constant times eps: a multiple of the identity. That is
why the whole near-null space moves rigidly. So to first order
"primes zoomed by eps" = "true form on a slightly different window, plus
a constant times eps, plus a pole term". The constant has a sign, which
is why one direction is allowed and the other is not; and the tiny
eigenvalues measure how much zoom-out the true primes could tolerate.

## The smallest level as zoom slack (`zoom_slack.py`, small windows where doubles resolve it)

Zoom p -> p^(1+eps); critical eps where the lowest eigenvalue of a half
reaches zero, measured by bisection, against lowest eigenvalue / its own
first-order rate.

| c | half | lowest eigenvalue | rate per unit zoom-out | predicted critical zoom | measured | ratio |
|---|---|---|---|---|---|---|
| 2.6 | odd | 7.04e-4 | 0.960 | +7.333e-4 | +7.325e-4 | 0.999 |
| 3.3 | odd | 8.06e-7 | 0.949 | +8.490e-7 | +8.487e-7 | 1.000 |
| 4.2 | odd | 5.95e-11 | 0.934 | +6.372e-11 | +6.370e-11 | 1.000 |
| 2.6 | even | 3.31e-6 | -0.117 | -2.838e-5 | -2.817e-5 | 0.993 |
| 3.3 | even | 2.05e-9 | -0.183 | -1.116e-8 | -1.113e-8 | 0.997 |

- **Odd half: the smallest eigenvalue is the zoom-out the primes can
  take, almost one for one** (rate 0.93 to 0.96, close to the 1 expected
  from a logarithmic symbol).
- **Even half: the lowest state moves the other way** (it is the
  pole-like state): it limits zooming *in*, at lowest eigenvalue / 0.12
  to 0.18.
- So with the pole direction kept, the scale is pinned from both sides,
  and the allowed interval is [-even level / 0.15, +odd level / 0.95]:
  [-2.8e-5, +7.3e-4] at c = 2.6, [-1.1e-8, +8.5e-7] at c = 3.3. It closes
  as fast as the smallest eigenvalues fall.

## How sharply does a window fix the scale of the primes?

Zoom slack = smallest odd eigenvalue / rate. Smallest eigenvalues from
the stored ball-arithmetic runs (`ladder_track.json`, `ladder_fine.json`):

| c | smallest odd eigenvalue | smallest even eigenvalue | odd: factor per unit of c |
|---|---|---|---|
| 2.2 | 1.8e-2 | 2.3e-4 | |
| 3.2 | 2.2e-6 | 6.1e-9 | 8.2e3 |
| 4.2 | 5.8e-11 | 9.7e-14 | 3.8e4 |
| 5.2 | 1.1e-15 | 9.8e-19 | 5.1e4 |
| 6.2 | 1.3e-20 | 9.5e-24 | 8.6e4 |
| 7.2 | 1.8e-25 | 6.0e-29 | 7.6e4 |
| 8.2 | 1.5e-30 | 5.0e-34 | 1.2e5 |

The factor per unit of c climbs towards e^(4 pi) = 2.9e5, the decay rate
of the Slepian (prolate) deficiency that Connes-Consani and ACL tie the
tiny eigenvalues to. So:

**a window up to c fixes the overall scale of the primes to a relative
precision of about e^(-4 pi c) times a slowly varying factor: 1e-6 at
c = 3, 1e-15 at c = 5, 1e-30 at c = 8.** The play in the scale is the
uncertainty-principle leakage of the window.

Status: the equality slack = eigenvalue / 0.95 is measured only up to
c = 4.2. Beyond, the rate is not measured for the lowest state itself;
on the whole near-null space at c = 14.2 the rates lie between 1.4 and
1.8 (pole: 3.7), so the slack is the smallest eigenvalue within a factor
of about 4. The zoom-in side (even half) is measured to c = 3.3 only and
involves the pole state; no law claimed.

## First look at the pattern stiffness (arithmetic on existing numbers, 26 Nov)

Guess: moving prime p shifts its terms at k ln p by k times the
displacement, so the natural size of the derivative is
S_p = sum over powers p^k <= c of k ln p / p^(k/2), times the highest
frequency the near-null states carry, 2 pi c.

| c | p | S_p | 2 pi c S_p | measured k | ratio |
|---|---|---|---|---|---|
| 14.2 | 2 | 1.92 | 171 | 47 / 37 | 0.25 |
| 14.2 | 3 | 1.37 | 122 | 25 / 27 | 0.21 |
| 14.2 | 7 | 0.735 | 66 | 14 / 13 | 0.21 |
| 14.2 | 11 | 0.723 | 64 | 13 / 10 | 0.18 |
| 14.2 | 13 (newest) | 0.711 | 63 | 7.7 / 5.2 | 0.10 |
| 28.2 | 3 | 2.00 | 354 | 91 / 90 | 0.26 |
| 28.2 | 7 | 0.735 | 130 | 53 / 39 | 0.35 |
| 28.2 | 23 (newest) | 0.654 | 116 | 24 / 18 | 0.18 |
| 56.2 | 7 | 1.29 | 455 | 123 / 115 | 0.26 |

So k_p is about a quarter of 2 pi c S_p, within a factor 1.5, except for
a prime that has just entered. That accounts for both trends (small
primes have more powers in the window; everything scales with the band
limit 2 pi c). The "c^1.5" seen earlier was three points and is replaced
by this. Not derived: the factor of about 1/4, and its scatter.

## The stiffness explained (`stiff_profile.py`)

The term of a prime power at log-position a is -w A_a, with
A_a(f, g) = int f(u) g(u+a) + g(u) f(u+a) du over the stretch of the
window where both points are inside. Checked against the form's own
matrices: weight derivative = -w A_a to 1e-12; position derivative
= -w (dA_a/da + (1/a - 1/2) A_a) to 1e-5 on the near-null space V.

So, exactly:

- **weight tent of a prime = w x extreme eigenvalues of A_a on V**
- **position tent = w x extreme eigenvalues of dA_a/da on V** (one power
  in the window; the weight part is negligible)

Predicted against measured position tents: c = 14.2: 5: 19.8 / 16.5
(measured 19.8 / 16.5); 7: 14.8 / 13.6 (14.8 / 13.6); 11: 13.8 / 11.4
(13.8 / 11.4); 13: 9.8 / 7.4 (9.8 / 7.4). c = 28.2: 7: 52.9 / 40.3
(53 / 39); 23: 24.8 / 19.6 (24.8 / 19.6).

**Weights.** The extremes of A_a on V are, to three or four digits,
+-2 cos(pi / (n + 1)) with n = number of steps of length a that fit in
the window: 1.000 for a > L/2, 1.414 for L/3 < a < L/2, 1.618 for
L/4 < a < L/3, 1.732, ... These are the norms of a hop on a chain of n
sites (the truncated shift). So the near-null space contains the
extremal vectors of the shift, and a prime beyond sqrt(c) has weight tent
exactly ln p / sqrt p. Only in the last fifth of the window does the
value drop below 1 (0.8 to 0.65 at a/L = 0.9 to 0.95): the newest primes
are not yet fully felt, the delay of `string.md`.

**Positions.** max |dA_a/da| on V, in units of 2 pi c:

| a / L | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 0.95 |
|---|---|---|---|---|---|---|---|---|---|---|
| c = 14.2 | 0.53 | 0.56 | 0.65 | 0.62 | 0.43 | 0.34 | 0.28 | 0.29 | 0.25 | 0.23 |
| c = 28.2 | 0.86 | 0.60 | 0.54 | 0.61 | 0.50 | 0.42 | 0.46 | 0.22 | 0.33 | 0.20 |

A falling, bumpy profile, similar at both windows: about 0.6 to 0.8
near the centre down to 0.2 to 0.3 at the far end. The "quarter" of the
previous section is its value over the outer half. So the stiffness of a
prime is w x 2 pi c x g(ln p / ln c), with g this profile; small primes
are stiffer because g is larger there and because their powers add.
The profile is a property of the near-null states alone (how much of
them, and how fast-varying, sits a apart), not of the primes.

## Chain model: tried, failed by the test set beforehand (`chain_model.py`, c = 14.2)

Goal stated first: ~25 sites along the window, each prime a link between
sites ln p apart, background on-site; success = at least 80% of each
prime's action on links of the right length and a mostly diagonal
background.

Sites = eigenvectors of the position operator compressed to the near-null
space (the most localised basis it allows).

| prime power | share on links of the right length | same in a random basis |
|---|---|---|
| 2 | 0.45 | 0.12 |
| 3 | 0.39 | 0.02 |
| 4 | 0.33 | 0.04 |
| 5 | 0.32 | 0.01 |
| 7 | 0.24 | 0.01 |
| 9 | 0.24 | 0.00 |
| 8, 11, 13 | 0.00 | 0.00 |

Background: 53% on the diagonal, 16% on nearest neighbours.

**Failed.** Far better than chance, far short of 80%. Reasons visible in
the output:

- The sites are not localised: each site's wave is 0.27 to 0.75 wide,
  against a spacing of 0.04 to 0.19. Twenty-five smooth functions cannot
  be turned into twenty-five points.
- The sites only reach +-1.08 in a window of +-1.33, so nothing is 2.08
  or more apart: the largest prime powers (8, 11, 13) link the two very
  ends of the window, beyond the outermost sites.
- The "hop on a chain of n sites" found for the weights is true prime by
  prime (slices of the window a apart), but each prime has its own
  chain, with a different step. There is no single set of sites that
  serves all of them, because the steps ln 2, ln 3, ln 5... are
  incommensurable.

So no small hopping model in position. If there is a small model it is
not labelled by place.

## Chain model, second attempt: sites labelled by frequency (`pitch_model.py`, c = 14.2)

Test stated first: at least 80% of each prime power's action on-site, and
on-site values within 10% of 2 (1 - a/L) cos(nu a).

| | even half | odd half |
|---|---|---|
| background on-site | 0.87 | 0.95 |
| 2 | 0.80 | 0.68 |
| 3 | 0.76 | 0.67 |
| 4 | 0.43 | 0.72 |
| 5 | 0.60 | 0.61 |
| 7 | 0.57 | 0.38 |
| 8, 9 | 0.51, 0.44 | 0.34, 0.33 |
| 11 | 0.32 | 0.28 |
| 13 | 0.15 | 0.13 |

On-site values miss the cosine by 0.3 to 0.5 on a range of +-1.

**Failed as stated, but half of it works.**

- The background is on-site in frequency (87 to 95%), with values rising
  like a logarithm from negative (-1.7, -0.5 at the lowest sites) to +2.3.
  That is the archimedean symbol, seen directly.
- Small primes are mostly on-site (2 and 3: 67 to 80%). Large primes are
  not (13: 15%). Their cosine has period 2 pi / ln p, about 2.4 to 3,
  which is the spacing of the sites: it cannot be resolved as a value
  per site.
- The odd sites at 2.4, 4.8, 7.2, 9.6, 12.0 are single Fourier modes, all
  below the first zero 14.13. Sites reach 60 only, of 89. They do not sit
  mid-gap between zeros.

Together with the position attempt: small primes are frequency-like
(a potential), the primes near the size of the window are place-like
(they join its two ends), and no single labelling makes all of them
simple. The near-null space is a phase-space object: about 25 cells,
neither places nor pitches.

## The "jammed packing" picture: tested, fails (`jam_test.py`, c = 14.2)

Speculation: the primes as a jammed pile at the edge of coming loose
(tents, no allowed neighbours, one soft breathing mode, a crowd of
nearly-free modes). Two standard signatures of jamming:

**(a) One contact removed, one mode goes soft.** Not so:

| left out | negative eigenvalues |
|---|---|
| one prime power (2, 3, 4, 5, 7, 8, 9, 11) | 8 to 13 |
| 13 alone (the newest) | 5 |
| 2 with 4 and 8 | 12 |
| 5 and 7 | 8 |
| 7, 11 and 13 | 11 |
| every prime | 4 |

Removing one prime power sends about half the near-null family negative,
the same as moving one. The count does not add up contact by contact,
and with no primes at all only 4 directions are negative (to -2.1). So
the primes are not local contacts: each one touches the whole family.

**(b) Soft modes with a flat density in frequency** (number below
stiffness k growing like sqrt k). Not so: the tiny levels at c = 8.2 are
1e-33.3, 1e-26.7, 1e-20.6, 1e-15.3, 1e-10.6, 1e-6.5, 1e-2.9 (even), in
steps of 4 to 7 decades. A geometric ladder, uniform in the logarithm.

**Verdict: not a jammed solid.** The tent and the rigidity are what any
system with many levels at the same energy does under a small push; they
do not single out jamming. A geometric ladder is the signature of
Boltzmann factors of evenly spaced modes (or of leakage through an ever
thicker barrier), which is the horizon / thermal reading already found
for the top of the ladder, and Eisler-Peschel's entanglement spectrum.
New fact from the test: the background alone has only 4 negative
directions at c = 14.2.

## The clock's ticks (`clock_ticks.py`, c = 8.2: primes 2, 3, 5, 7 only; N = 50, 60 digits)

Connes's recipe: lowest even state xi of the form, its transform
xi^(s) = int xi(u) cos(su) du, its real zeros ("ticks": the spectrum of
the scaling operator on the window with the boundary rule built from xi).

**True primes: 13 ticks below 60, and zeta has 13 zeros there.**

| zero | 14.13 | 21.02 | 25.01 | 30.42 | 32.94 | 37.59 | 40.92 | 43.33 | 48.01 | 49.77 | 52.97 | 56.45 | 59.35 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| tick - zero | 3e-30 | 4e-27 | 2e-25 | 1e-22 | 2e-21 | 1e-19 | 7e-18 | 6e-17 | 8e-15 | 6e-14 | 4e-13 | 1e-11 | 4e-10 |

Nothing below 14.13. (This reproduces CCM / Connes-van Suijlekom; the
zeros from `mpmath.zetazero`.)

**Primes moved:**

| move | real ticks below 60 | the genuine ones | extra ticks |
|---|---|---|---|
| none | 13 | exact to 1e-30 ... 1e-10 | none |
| zoom out 0.1% | 17 | at zero / 1.001, to 1e-5 ... 1e-2 | 3 below the first zero, 1 near 38 |
| zoom out 1% | 17 | at zero / 1.01, to 1e-3 ... 1e-2 | 4.22, 7.56, 11.67 and 21.94 |
| zoom in 0.1% | 19 | at zero / 0.999, to 1e-4 ... 1e-2 | 3 below the first zero, 3 higher |
| 7 -> 7.007 | 19 | shifted by 1e-4 to 3e-2 | 0.54, 3.44, 6.76, 10.14, 16.61, 24.00 |
| 7 -> 7.07 | 19 | shifted by 2e-3 to 0.7 | 1.06, 3.77, 6.93, 10.23, 16.67, 23.88 |
| 2 -> 2.02 | 19 | shifted by 0.03 to 0.5 | 1.53, 4.60, 7.50, 10.54, 19.79, 27.54 |

- **Under a zoom the ticks follow the tape speed**: they sit at
  zero / (1 + eps). The clock really is the zoom.
- **Any move fills the silence below the first zero with extra ticks,
  about 3.0 to 3.4 apart.** The empty box's own notes are 2 pi / L = 2.99
  apart. So with the primes even slightly wrong, the clock ticks like the
  empty box at low frequency; with the true primes those notes are
  cancelled exactly and the first tick is at 14.13.
- The extra ticks appear at once (the same set at 0.1% and at 1%), not
  gradually: the lowest state has changed identity (its eigenvalue is
  -0.005, not 6e-34).
- All ticks found are real in every case; the count goes up, never down.
  Complex ticks were not searched for.

### Following the true vacuum instead of the lowest state (`FOLLOW=1 clock_ticks.py`, c = 8.2)

For each move: the eigenstate of the moved form with the largest overlap
with the true lowest state, and that state's ticks.

| move | position of the prime in the window (ln p / L) | overlap with the true vacuum (next best) | its eigenvalue, rank from bottom | ticks below 60 | error of the first zero |
|---|---|---|---|---|---|
| 7 by 0.1% | 0.92 | 0.99991 (0.014) | 1e-24, 4th | 13, no extras | 4e-21 |
| 7 by 1% | 0.92 | 0.99992 (0.012) | -4e-23, 4th | 13 | 1e-19 |
| 5 by 0.1% | 0.76 | 0.9972 (0.07) | -1e-10, 4th | 13 | 1e-7 |
| 5 by 1% | 0.76 | 0.9967 (0.08) | -8e-9, 4th | 13 | 7e-6 |
| 3 by 0.1% | 0.52 | 0.79 (0.55) | 8e-4, 5th | 19 | lost |
| 3 by 1% | 0.52 | 0.80 (0.55) | 8e-3, 5th | 17 | lost |
| 2 by 0.1% | 0.33 | 0.92 (0.33) | 8e-4, 4th | 17 | lost |
| 2 by 1% | 0.33 | 0.80 (0.56) | 9e-3, 5th | 18 | lost |
| zoom out / in 1% | all | 0.80 / 0.82 (0.5) | -6e-3 / +6e-3 | 18 | scaled, with extras |

- **Moving an outer prime (7, 5) leaves the vacuum intact.** It is still
  there, almost unchanged, just no longer the lowest state; its ticks are
  the 13 zeros with no extras. With 7 at 7.07 the first zero is still
  right to 19 digits. The breakdown seen before was in the choice of
  state, not in the state.
- **Moving an inner prime (2, 3), or zooming, destroys it.** No single
  state resembles the vacuum; it is spread over several, and the extra
  ticks are real.
- The dividing line is position in the window: the vacuum is small near
  the ends, and a prime near the size of the window only compares the two
  ends. This is the "flat, edge-blind" state of `ladder.md` and the edge
  response of `string.md` again.
- Consequence: the low zeros are fixed to extreme precision by the small
  primes and the background alone. The largest primes in the window
  contribute almost nothing to them (first zero: about 1e-17 per unit
  displacement of 7, 7e-4 per unit of 5).

## The map, version 1 (`sens_map.py`, c = 8.2, N = 50, 90 digits, delta = 1e-36)

Cell = log10 of |d(tick) / d(delta)| when prime p is moved to p (1 + delta),
in the linear regime (delta far below the gap 2.8e-27 between the vacuum
and the next state, so the vacuum stays the lowest state; its overlap with
the true one differs from 1 by 1e-19 to 1e-49).

| zero | 2 | 3 | 5 | 7 |
|---|---|---|---|---|
| 14.13 | +2.7 | +1.3 | -4.5 | -13.8 |
| 21.02 | +5.4 | +4.1 | -1.5 | -10.7 |
| 25.01 | +7.0 | +5.8 | +0.3 | -8.9 |
| 30.42 | +9.5 | +8.4 | +2.9 | -6.3 |
| 32.94 | +10.6 | +9.5 | +4.0 | -5.1 |
| 37.59 | +12.4 | +11.2 | +5.9 | -3.3 |
| 40.92 | +14.0 | +12.9 | +7.6 | -1.5 |
| 43.33 | +14.9 | +13.8 | +8.4 | -0.6 |
| 48.01 | +16.9 | +15.8 | +10.5 | +1.5 |
| 49.77 | +17.7 | +16.6 | +11.4 | +2.3 |
| 52.97 | +18.5 | +17.4 | +12.2 | +3.1 |
| 56.45 | +19.8 | +18.7 | +13.5 | +4.5 |
| 59.35 | +21.3 | +20.2 | +15.0 | +6.0 |

**Not a staircase: a product.** Column differences are constant down the
table (2 minus 3: 1.1 to 1.4; 3 minus 5: 5.2 to 5.8; 5 minus 7: 9.0 to
9.3). So cell = (touchiness of the zero) x (strength of the prime):

- strength of the prime, relative to 2: 3: 10^-1.1, 5: 10^-6.5, 7: 10^-15.6.
  It collapses towards the edge of the window.
- touchiness of the zero: grows by 18.6 decades from the 1st to the 13th
  zero, about 0.4 decades per unit of height (e^(pi gamma / 4) would give
  0.34: the decay rate of the xi function).

Why a product: in the linear regime every prime disturbs the vacuum the
same way, by mixing in its one nearest neighbour (eigenvalue 2.8e-27).
The prime decides how much is mixed in; the zero decides how much that
admixture moves it. One direction of disturbance, for all primes.

Reading: to keep the 13th zero right to 1e-10 one needs 2 placed to
1e-31, 5 to 1e-25, 7 to 1e-16. These tolerances are of the order of the
tiny eigenvalues: the same slack as for the zoom.

Caution: this is the linear regime, valid for displacements below about
1e-27. For large displacements of an outer prime the response is far
smaller than linear (7 moved by 1e-2: first zero off by 1e-19, linear
would say 1e-16), because the vacuum stops being the lowest state and
decouples. The two regimes answer different questions.

## The map, version 2: three windows (`sens_map.py`, delta = 1e-70, 150 digits)

c = 14.2 at N = 70 failed numerically (lowest eigenvalues came out
-1.7e-55 and -1.9e-56: the quadrature is not accurate enough at that
depth; needs HP_DIGITS raised). Discarded. Windows 5.2, 8.2, 11.2:

**The product structure holds at every window.** Column differences are
constant down each table to a few tenths.

Strength of each prime relative to 2 (log10), with its position ln p / L:

| prime | c = 5.2 | c = 8.2 | c = 11.2 |
|---|---|---|---|
| 3 | -1.5 (0.67) | -1.1 (0.52) | -1.0 (0.46) |
| 5 | -10.5 (0.98) | -6.5 (0.77) | -5.8 (0.67) |
| 7 | | -15.6 (0.93) | -12.5 (0.81) |
| 11 | | | -39 (0.99) |

Sensitivity of the first zero to prime 2 (log10 per unit displacement):
+2.5, +2.7, +2.8. Of the 13th tick: +14.6, +21.3, +23.6.

- **A prime that has just entered is fantastically weak** (5 at c = 5.2:
  1e-10.5; 11 at c = 11.2: 1e-39), **and gains strength as the window
  moves past it** (5: -10.5, -6.5, -5.8; 7: -15.6, -12.5). The inner
  primes' strengths settle (3: -1.5, -1.1, -1.0).
- **The first zero's sensitivity to prime 2 is stable, about 300 to 600
  per unit**, although the gap above the vacuum falls from 8e-13 to
  3e-42 over these windows. So the huge mixing with the neighbouring
  state (amplitude about delta / gap) is almost exactly cancelled: the
  neighbour's transform also nearly vanishes at the zeros, being a
  near-null state itself. The product structure is consistent with "one
  admixture, for all primes", but the size of the cells is not set by
  the gap. Not understood beyond that.
- The touchiness of the higher ticks grows with the window (13th: 14.6,
  21.3, 23.6), as more of them become genuine zeros (2 pi c = 33, 52, 70).

A caution about meaning: for a true Euler product, moving one prime
would not move the zeros of zeta at all (it changes one factor). These
sensitivities belong to the recipe on a window: they measure how the
vacuum of the *moved* form, and hence its boundary rule, responds.

## The map, version 2b: window 14 fixed, and the fade-in of each prime

Window 14.2 fixed with `HP_DIGITS=80` (the default graded quadrature
stops at 2^-180 = 7e-55; the vacuum there is 2.8e-65, the next state
1.7e-57). Product structure holds there too.

Strength of each prime relative to 2, from the first zero's row (log10),
windows 3.3 to 14.2, all with `HP_PREC=600 HP_DIGITS=80`, delta = 1e-70:

| c | 3 | 5 | 7 | 11 | 13 |
|---|---|---|---|---|---|
| 3.3 | -3.3 | | | | |
| 4.2 | -2.2 | | | | |
| 5.2 | -1.9 | -11.9 | | | |
| 6.2 | -1.7 | -8.8 | | | |
| 7.2 | -1.6 | -7.8 | -20.9 | | |
| 8.2 | -1.4 | -7.2 | -16.5 | | |
| 9.2 | -1.4 | -6.9 | -14.9 | | |
| 11.2 | -1.3 | -6.4 | -13.2 | -39.9 | |
| 14.2 | -1.3 | -6.1 | -12.2 | -28.8 | -42.1 |

First zero's sensitivity to prime 2 itself: +2.1, +2.3, +2.5, +2.6, +2.7,
+2.7, +2.8, +2.8, +2.9 (slowly rising, about 125 to 800 per unit).

**Is the fade-in one curve in c / p? No.**

- Each prime settles to its own level, and the levels fall steeply with
  the prime: about -1.3 (3), -6 (5), -12 (7), and 11 is still falling in
  at -29. That is 2.4 to 3 decades per unit of p. (e^(2 pi) is 2.73
  decades. Possibly the same leak rate as the tiny levels, per unit of p
  instead of per unit of c; three points, not a claim.)
- On entry a prime is weaker again by about the same number of decades
  as its settled level: 5 enters at -11.9 and settles near -6; 7 enters
  at -20.9 and settles near -12; 3 enters at -3.3 and settles at -1.3.
  So roughly log(strength) = (settled level of p) x g(c / p), with g about
  1.8 to 2.5 just after entry, 1.4 at c / p = 1.2, 1.2 at 1.5, 1 beyond 3.
  This holds to 10 to 20% for 5 and 7 and worse for 3.

So: **every prime fades in over about a doubling of the window, and
bigger primes settle at exponentially smaller strengths.** In this
recipe the low zeros are made by 2 and 3; 5 contributes a millionth, 7 a
trillionth.

## The map, version 3: larger windows by perturbation theory (`sens_fast.py`)

One build of the form, its inverse, the vacuum by inverse iteration, and
the exact derivative of the form from the shift operators in closed
form. Reproduces `sens_map.py` at c = 14.2 to the digits printed
(+2.87, +1.61, -3.16, -9.27, -25.92, -39.17) in 7 s instead of 65.

| c | vacuum eigenvalue | first tick minus first zero | 2 | 3 | 5 | 7 | 11 | 13 | 17 | 19 | 23 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 14.2 | 2.8e-65 | 2.1e-61 | +2.87 | +1.61 | -3.16 | -9.27 | -25.92 | -39.17 | | | |
| 18.2 | 3.6e-85 | 3.0e-81 | +2.92 | +1.71 | -2.83 | -8.48 | -22.35 | -30.97 | -56.54 | | |
| 22.2 | 7.6e-104 | 6.9e-100 | +2.95 | +1.78 | -2.63 | -8.01 | -20.66 | -27.99 | -45.69 | -57.22 | |
| 28.2 | 4.3e-130 | 4.3e-126 | +2.98 | +1.84 | -2.44 | -7.58 | -19.27 | -25.74 | -40.14 | -48.29 | -67.68 |

(log10 of |d(first tick)/d(delta)| for each prime.)

**1. The test of the tentative link: it fails.** There is no settled
level and no constant fall per unit of p. Every prime keeps gaining
slowly (5: -6.0, -5.8, -5.6, -5.4 relative to 2), and at c = 28.2 the fall
per unit of p rises from 2.1 (3 to 5) to 4.9 (19 to 23). "2.73 decades
per unit of p" was three points.

**2. What does hold: one curve, in the right variables.** Write
t = ln p / ln c (where the prime sits in the window) and measure the
strength relative to 2 in units of the vacuum's own depth:
h = log10(strength_p / strength_2) / log10(vacuum eigenvalue).

| t | 0.33 | 0.41 | 0.48 | 0.52 | 0.58 | 0.63 | 0.67 | 0.72 | 0.77 | 0.83 | 0.88 | 0.91 | 0.95 | 0.98 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| h | 0.009 | 0.020 | 0.042 | 0.054 | 0.082 | 0.106 | 0.135 | 0.172 | 0.229 | 0.300 | 0.40 | 0.47 | 0.58 | 0.70 |

Points from all four windows fall on it. Same t from different windows:
t = 0.827: 0.300 (13 at c = 22.2) and 0.299 (11 at c = 18.2). t = 0.88:
0.401 (13 at 18.2) and 0.396 (19 at 28.2). t = 0.77: 0.229 (11 at 22.2)
and 0.222 (13 at 28.2).

So **strength of a prime = (vacuum eigenvalue)^h(t)**, relative to 2: a
prime at the very edge is as weak as the vacuum is deep, and one in the
inner third has nearly full strength. The fade-in is a property of the
position on the log ruler, the same for every prime and every window.

**3. A by-product: the tiny level is the accuracy of the zeros.** The
first tick misses the first zero by 7.3e3, 8.3e3, 9.1e3, 9.9e3 times the
vacuum eigenvalue at the four windows. With primes up to 23 the recipe
gives the first zero to 126 digits.

The first zero's sensitivity to prime 2 itself creeps up: +2.87, +2.92,
+2.95, +2.98.

## The formula for the fade-in: a prime reads the vacuum at sqrt(p)

**Result.** For a prime with one power in the window,

    strength_p  ~  (ln p / sqrt p) x p^nu x exp( -4 pi c (1 - sqrt(1 - p/c)) )

with nu about 8 to 11 (fitted). The exponential alone removes all but a
smooth remainder: over 21 points (c = 14.2 to 28.2, p = 5 to 23) the
sensitivities span 65 decades (-3 to -68); after dividing by the weight
and the exponential the remainder runs from 12 to 20 decades and is, to
0.1 decade at the larger windows, a function of p alone
(5: 12.03; 7: 13.01; 11: 14.58; 13: 15.32; 17: 16.93). With a power of p
the fit is good to 0.5 decade rms.

**Where it comes from (a sketch, not a derivation).**

1. In the additive variable y (window [-lambda, lambda], c = lambda^2)
   the vacuum is a prolate-type wave. Its WKB envelope, from the prolate
   equation with small eigenvalue, is
   phi(y) ~ exp( -2 pi c (1 - sqrt(1 - y^2/c)) ):
   a bell exp(-pi y^2) in the middle, falling to exp(-2 pi c) at the end.
2. The prime p compares the state at x with the state at p x. By the
   symmetry x -> 1/x that is phi(y) phi(p / y), with y between
   p / lambda and lambda.
3. The product is largest at y = sqrt(p), where it is
   exp( -4 pi c (1 - sqrt(1 - p/c)) ).

So **the influence of a prime is the intensity of the vacuum at
sqrt(p)**. Limits: for p far inside, exp(-2 pi p) (the "2.73 decades per
unit of p" seen at first, which is right only there); for p at the edge,
exp(-4 pi c), the depth of the vacuum itself (h -> 1). The collapse on
h(ln p / ln c) of the previous section is this formula, approximately;
the exact scaling variable is p / c.

Not derived: the power nu (it should come from the prolate index of the
vacuum, n + 1/2 per factor, and from the WKB prefactors; the fitted 8 to
11 is larger than n = 6 would give), and primes with several powers in
the window.

## Pinning down the power of p (`fade_fit.py`)

Instead of fitting a power, use the WKB wave of the prolate equation with
its eigenvalue kept, chi ~ 2 pi c (2n + 1), prefactor included:

    psi_n(y) ~ [ (4 pi^2 y^2 - 2 pi (2n+1)) (1 - y^2/c) ]^(-1/4)
               x exp( - int sqrt( (4 pi^2 x^2 - 2 pi (2n+1)) / (1 - x^2/c) ) dx )

and model strength_p = (ln p / sqrt p) psi_a(sqrt p) psi_b(sqrt p) x one
constant per window. Scan of the index pairs (a, b), 21 points spanning
65 decades:

| indices | rms of the remainder (decades) |
|---|---|
| 6, 6 (best; any pair with a + b = 12 is as good) | 0.50 |
| 2, 6 | 0.80 |
| 2, 2 | 1.32 |
| 0, 0 (the plain bell, no power) | 1.86 |

**The power is that of the sixth prolate wave.** Only the sum of the two
indices is determined (12), but 6 is the lowest index of the even ladder
found in `shape_brainstorm.md` (n = 6, 10, 14, ...), so the natural
reading is: the vacuum is the n = 6 member, and a prime reads its
intensity at sqrt(p). In the bulk that is about p^6 exp(-2 pi p); the
freely fitted power (8 to 11) was larger because near its turning point
the wave is not a pure power.

Remainders are within 0.7 decade except 23 at c = 28.2 (+1.6, close to
the edge), and they still drift with the window (constants 2.9, 3.2,
3.6, 4.4), so this is a good leading description, not an exact formula.

## Extending the formula to 2 and 3: one response curve for a mirror at any distance (`response_curve.py`, c = 28.2)

For a shift a = ln q, q any real number in the window:
TA(q) = change of the first tick when a unit mirror at q is strengthened
(perturbation -A_a), T(q) = when it is moved (perturbation -A'_a). Then
for any prime, exactly,
d tick / d delta = sum_k p^(-k/2) [ k ln p T(p^k) + (1 - (k/2) ln p) TA(p^k) ].

| q | 1.1 | 1.3 | 1.5 | 1.7 | 1.9 | 2 | 2.2 | 2.5 | 3 | 3.5 | 4 | 5 | 6 | 7 | 8 | 9 | 11 | 13 | 16 | 17 | 19 | 23 | 25 | 27 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| log10 of -TA | 2.08 | 2.82 | 2.96 | 2.90 | 2.73 | 2.62 | 2.34 | 1.86 | 0.90 | -0.16 | -1.31 | -3.76 | -6.37 | -9.10 | -11.93 | -14.87 | -21.03 | -27.60 | -38.32 | -42.17 | -50.40 | -69.97 | -82.61 | -100.34 |
| T (sign, log10) | -3.38 | -3.49 | -2.48 | +3.29 | +3.40 | +3.37 | +3.22 | +2.86 | +2.05 | +1.09 | +0.03 | -2.29* | -4.80* | -7.44* | -10.21* | -13.08* | -19.13* | -25.59* | -36.17* | -39.98* | -48.12* | -67.49* | -80.00* | -97.48* |

(* positive, the number is log10 T.)

1. **One smooth curve, primes or not.** The prime powers 4, 8, 9, 16, 25,
   27 lie on the same curve as the primes. A mirror at distance q is read
   the same way whatever q is.
2. **For 2 and 3 the first power is everything.** Prime 2: 956 from q = 2,
   0.73 from 4, 5e-11 from 8, 5e-37 from 16. Prime 3: 69 from 3, 6e-14
   from 9. So the extension is immediate: strength_2 = (ln 2 / sqrt 2) T(2)
   + small, and the same formula covers every prime.
3. **TA is an intensity and T is its slope.** TA never changes sign and
   peaks at q about 1.5; T changes sign there (between 1.5 and 1.7) and
   T / TA equals the logarithmic slope of TA (q = 6: 37 against 36.5;
   q = 11 to 13: 79 to 102 against 91). So the tick depends on the primes
   through one sum:

       first tick = constant - K sum over prime powers q of (ln p / sqrt q) I(q)

   to first order, with I(q) = -TA(q)/K the kernel.
4. **The kernel is the sixth Hermite wave's intensity at sqrt(q), with
   the edge factor.** I(q) ~ H_6(sqrt(2 pi q))^2 exp(-4 pi c (1 - sqrt(1 -
   q/c))): the offset between this and log10(-TA) is 1.88, 1.44, 1.48,
   1.55, 1.83 at q = 1.1, 1.5, 2, 3, 5 (constant to +-0.25 while TA varies
   by 6.7 decades), and its outermost maximum is at q about 1.5, where TA
   peaks and T changes sign. The plateau for q below about 2.1 is the
   wave's oscillating region (turning point sqrt(13 / 2 pi), q = 2.07);
   2 sits just inside it, 3 just outside. That is why 2 and 3 carry
   nearly everything.

So prime 2 is special for a simple reason: sqrt(2) is the only prime
root that falls inside the bright part of the vacuum.

### Checks of the kernel: other zeros, another window (`TICKS=... response_curve.py`, c = 14.2 and 28.2)

**Other zeros (c = 28.2, ticks 1 to 6).** log10 of TA_n / TA_1 across q:

| q | 1.3 | 2 | 3 | 5 | 7 | 11 | 17 | 23 |
|---|---|---|---|---|---|---|---|---|
| tick 2 | 3.04 | 3.13 | 3.26 | 3.38 | 3.42 | 3.45 | 3.46 | 3.46 |
| tick 4 | 8.00 | 8.12 | 8.31 | 8.54 | 8.65 | 8.72 | 8.75 | 8.76 |
| tick 6 | 11.67 | 11.79 | 12.01 | 12.30 | 12.44 | 12.56 | 12.60 | 12.62 |

Beyond q about 7 the ratio is constant to a few hundredths while TA falls
by 60 decades: there the kernel is the same for every zero. In the bright
region (q below 5) it tilts by 0.3 to 0.6 decade, more for higher zeros.
So "one kernel for all zeros" is exact in the dark part and approximate
in the bright part, where 2 and 3 sit.

**Another window.** log10(-TA) for the first tick, c = 14.2 minus
c = 28.2, against the difference of the edge factors
4 pi c (1 - sqrt(1 - q/c)) / ln 10:

| q | 2 | 3 | 5 | 7 | 13 |
|---|---|---|---|---|---|
| measured difference | -0.15 | -0.26 | -0.78 | -1.79 | -13.85 |
| edge factors predict | -0.10 | -0.25 | -0.81 | -1.85 | -14.06 |

The window enters only through the edge factor, to 0.05 to 0.2 decade
over 14 decades. Same peak (q about 1.5) and same sign change of T at
both windows.

### The small tilt (`DECOMP=J response_curve.py`, c = 28.2)

Decompose the response over the eigenstates above the vacuum:
TA_n(q) = sum_j c_j(q) m_jn, with c_j(q) = <v_j | A_a | v0> and
m_jn = v_j^(t_n) / (lambda_j v0^'(t_n)). Eigenvalues: 4.3e-130 (vacuum),
9.8e-122, 4.7e-114, 7.5e-107, 2.8e-100, ...

**Bright region (q up to about 2): two states, and that is the tilt.**

| q | tick 1: share of j = 1, j = 2 | tick 6: share of j = 1, j = 2 |
|---|---|---|
| 1.5 | +1.12, -0.12 | +1.01, -0.01 |
| 2 | +1.43, -0.54 | +1.03, -0.04 |

The second state enters with the opposite sign, and couples to the first
zero ten times more strongly than to the sixth (m_2n / m_1n = 10^-0.6
for tick 1, 10^-1.6 for tick 6). Predicted change of TA_6 / TA_1 from
q = 1.5 to 2: log10((1.43/1.03) / (1.12/1.01)) = +0.10. Measured: +0.09.
So the tilt is the second near-null state: it matters for the lowest zero
and fades for higher ones.

**Dark region (q beyond about 3): the state-by-state picture breaks
down.** The coefficients c_j(q) / c_1(q) grow steeply with j and q (at
q = 23: 1e3, 1e6, 1e9, 1e11, ...), the terms alternate in sign, and the
partial sum over the first 12 states is 10^-0.4 at q = 7 where the true
value is 10^-9.1. Many states cancel to tens of digits. There the
response is not "one admixture"; it is a collective cancellation whose
result happens to be the smooth kernel, the same for every zero.

This corrects the earlier explanation of the product structure ("every
prime mixes in the one nearest state"): true for q below about 2, false
beyond. Why the collective sum is nevertheless a clean product in the
dark region is not understood.

### Why the collective sum is a clean product: a mechanism, and a test (`OSC=1 response_curve.py`)

**Argument (a sketch).**

1. The form is a sum over zeros: Q = 2 sum_gamma |e_gamma><e_gamma|,
   e_gamma(u) = cos(gamma u) on the window. The first-order shift of tick
   n under a perturbation X is <G_n, X v0>, where G_n = Q^-1 e_(t_n) /
   v0^'(t_n) is the form's response to a source at that zero.
2. Q G = e_(gamma_n) asks for a function whose transform vanishes at
   every zero except gamma_n. The vacuum's transform vanishes at all of
   them, so the natural solution is the vacuum with that zero divided
   out: G_n^(s) ~ v0^(s) / (s^2 - gamma_n^2), i.e.
   (-d^2/du^2 - gamma_n^2) G_n ~ v0. It stays inside the window exactly
   because v0^(gamma_n) = 0.
3. For a mirror at lag a (X = -A_a) this gives
   TA_n(a) ~ (1 / gamma_n) int_a^L R(y) sin(gamma_n (y - a)) dy,
   R(y) = int v0(u) v0(u + y) du the vacuum's self-overlap; equivalently
   TA_n'' + gamma_n^2 TA_n ~ R(a).
4. In the dark region R falls so steeply beyond a that sin(gamma_n (y -
   a)) is linear over the part that counts; the 1 / gamma_n cancels and
   the zero drops out: TA_n(a) ~ int_a^L (y - a) R(y) dy for every n. In
   the bright region R varies slowly, the sine matters, and the tilt
   appears, larger for higher zeros. The crossover is where the decay
   length of R equals 1 / gamma_n: q about 3 for the first zero, about 6
   for the sixth, as observed.

So the cancellation over many states is orderly because the states were
the wrong variables: the response function of a zero is one object, the
vacuum with that zero removed.

**Test.** (TA_n'' + t_n^2 TA_n) / R(a), c = 28.2, ticks 1 to 4:

| q | log10 R | tick 1 | tick 2 / tick 1 | tick 3 / tick 1 | tick 4 / tick 1 |
|---|---|---|---|---|---|
| 1.3 | -0.16 | 1.7e5 | 2604 | 3.32e5 | 5.03e8 |
| 2 | -1.15 | 1.3e6 | 2860 | 3.71e5 | 5.71e8 |
| 3 | -3.16 | 4.4e6 | 2889 | 3.76e5 | 5.79e8 |
| 5 | -7.98 | 1.7e7 | 2897 | 3.77e5 | 5.81e8 |
| 9 | -19.15 | 7.5e7 | 2900 | 3.78e5 | 5.82e8 |
| 23 | -74.26 | 1.8e9 | 2900 | 3.78e5 | 5.82e8 |

- **Putting the zero back removes most of the tilt.** Between q = 1.3 and
  the dark region the raw ratio TA_2 / TA_1 changes by a factor 2.6;
  after applying d^2/da^2 + t_n^2 it changes by 11% (tick 4: factor 5.6
  raw, 16% after). The oscillator structure is right to that accuracy.
- **The source follows the vacuum's self-overlap in exponential order
  but not exactly.** R falls by 74 decades; the quotient drifts by 4
  (about q^2.6, more near the edge). So step 2 holds up to a slowly
  varying factor. Likely reasons, not checked: the zeros above 2 pi c,
  which the vacuum's transform cannot all vanish at, and the pole terms
  of the form.

### The loose thread resolved: the drift was the vacuum's level (`KCHK=1 response_curve.py`, c = 28.2)

The slowly varying factor came from a term the sketch left out. From the
eigen-equation of the perturbed form, (Q + eps X) v = lambda v with
Q = 2 sum |e_gamma><e_gamma|, pairing with d_n (the vacuum with zero n
divided out) gives two terms:

    d(tick_n) = kappa_n [ <d_n, v0> <v0, X v0>  -  <d_n, X v0> ]
                          level term               shape term

The level term is K_n times the first-order shift of the vacuum's own
eigenvalue, with K_n = (tick_n - zero_n) / lambda_0 (the same pairing
applied to the unperturbed equation gives tick_n = zero_n + K_n lambda_0).
The shape term is the oscillator solution of the previous section. Since
the level term is proportional to R(a) itself, applying d^2/da^2 + t_n^2
to it brings in R''(a), which grows like (slope)^2: that was the "q^2.6"
drift (slope^2 ratio between q = 3 and 23: 450; drift ratio: 404).

**Test.** K_n from the unperturbed ticks: 9923.9, 2.885e7, 3.764e9,
5.818e12. Then TA_n / (K_n x d(level)), d(level) = -<v0|A_a|v0>:

| q | 1.5 | 2 | 3 | 5 | 7 | 9 | 13 | 17 | 23 | 27 |
|---|---|---|---|---|---|---|---|---|---|---|
| tick 1 | 0.111 | 0.298 | 0.589 | 0.837 | 0.919 | 0.954 | 0.981 | 0.991 | 0.996 | 0.998 |
| tick 2 | 0.044 | 0.139 | 0.367 | 0.686 | 0.831 | 0.900 | 0.957 | 0.979 | 0.991 | 0.995 |
| tick 4 | 0.020 | 0.066 | 0.204 | 0.498 | 0.695 | 0.808 | 0.913 | 0.955 | 0.982 | 0.990 |

**So: tick_n = zero_n + K_n x (vacuum level), and a far mirror does
nothing but change the level.** The response to a mirror at lag a is
K_n times -2 R(a), the vacuum's self-overlap, to 0.2% at the edge of the
window, over a range where R falls by 104 decades. The deficit from 1 is
the shape term, of relative size about (t_n x decay length of R)^2
(q = 13: (14.13 / 100)^2 = 0.020 against 0.019 for tick 1; (30.4 / 100)^2
= 0.092 against 0.087 for tick 4). For near mirrors the two terms nearly
cancel (0.11 at q = 1.5, 0.30 at q = 2 for the first zero): that
cancellation, different for each zero, is the tilt.

This replaces "the kernel is the sixth wave's intensity at sqrt(q)" by
the exact statement behind it: **the kernel is the vacuum's self-overlap
at lag ln q** (of which the wave formula is the WKB estimate), and the
"touchiness" of a zero is K_n, its tick's offset per unit of vacuum level.

## Near primes: the bright region is Riemann's own function (`near_primes.py`)

As the window grows the vacuum's transform tends to Xi(s) = xi(1/2 + is)
(CCM), so the vacuum itself tends to the classical function

    Phi(u) = 2 sum_n (2 pi^2 n^4 e^(9u/2) - 3 pi n^2 e^(5u/2)) exp(-pi n^2 e^(2u)),

whose cosine transform is Xi. In the variable y = e^u it is a bell
y^(9/2) exp(-pi y^2): the "bell with a power" of the fade-in formula,
without the window. Then, with no window and no fitted number:

- self-overlap R(a) = int Phi(u) Phi(u + a) du / int Phi^2 (transform: Xi^2);
- the same with zero gamma_n divided out, P_n(a) = (1 / gamma_n)
  int_a^inf R(y) sin(gamma_n (y - a)) dy (transform: Xi^2 / (s^2 - gamma_n^2));
- response of zero n to a mirror at lag a, in units of K_n times the
  level shift: rho_n(a) = 1 - (P_n(a) / P_n(0)) / R(a).

| q | rho from Phi alone: zeros 1, 2, 4 | measured at c = 28.2 |
|---|---|---|
| 1.5 | 0.106, 0.042, 0.019 | 0.111, 0.044, 0.020 |
| 2 | 0.283, 0.131, 0.063 | 0.298, 0.139, 0.066 |
| 3 | 0.564, 0.344, 0.189 | 0.589, 0.367, 0.204 |
| 5 | 0.813, 0.650, 0.459 | 0.837, 0.686, 0.498 |
| 7 | 0.901, 0.797, 0.646 | 0.919, 0.831, 0.695 |

Agreement to 3 to 8% everywhere, the difference being the finite window.
The self-overlap itself: log10 R(ln q) = -0.371, -1.121, -3.049 from Phi
against -0.38, -1.15, -3.16 measured.

**So the near-prime behaviour is not a property of our windows. It has a
limit, and the limit is written in terms of Xi alone:** the kernel is the
transform of Xi^2, and the cancellation between level and shape for a
given zero is the transform of Xi^2 with that zero's factor removed. For
the prime 2 and the first zero, 72% of the level effect is cancelled;
for the fourth zero, 94%.

What does not obviously have a limit: K_n itself (7350, 8340, 9110, 9924
for the first zero at c = 14.2, 18.2, 22.2, 28.2), multiplying a level
that goes to zero.

## The K multiplier (`k_multiplier.py`)

**Identity.** Pairing Q v0 = lambda_0 v0 with the vacuum-minus-one-zero
wave, with v0^ the cosine transform of the normalised vacuum:

    K_n = (tick_n - gamma_n) / lambda_0
        = (t_n / v0^'(t_n)^2) [ P_n(0) + 2 sum_{gamma != gamma_n} w_gamma / (gamma^2 - t_n^2) ]

where w_gamma = v0^(gamma)^2 / lambda_0 are the shares of the vacuum's
level carried by each zero (sum of 2 w = 1) and P_n(0) is about 1 / t_n^2
(t_n^2 P_n(0) = 1.061, 1.025, 1.017, 1.011 for zeros 1 to 4).

**Limit from Riemann's Phi alone** (first part, window to infinity):
K_n^Phi = gamma_n P_n(0) / Xi~'(gamma_n)^2, Xi~ the transform of Phi / |Phi|.

| zero | Xi~'(gamma) | K^Phi | K measured, c = 28.2 | ratio |
|---|---|---|---|---|
| 1 | -2.445e-3 | 1.256e4 | 9.924e3 | 0.790 |
| 2 | +3.139e-5 | 4.948e7 | 2.885e7 | 0.583 |
| 3 | -2.238e-6 | 8.124e9 | 3.764e9 | 0.463 |
| 4 | +4.265e-8 | 1.827e13 | 5.818e12 | 0.318 |

First zero, by window: K / K^Phi = 0.585, 0.665, 0.725, 0.790 at
c = 14.2, 18.2, 22.2, 28.2. Then (1 - ratio) x c = 5.9, 6.1, 6.1, 5.9.

**So K has a limit, the classical number gamma_n P_n(0) / Xi~'(gamma_n)^2
(essentially 1 / (gamma_n Xi'(gamma_n)^2)), approached like 1 - 6.0 / c
for the first zero.** The creep with the window was a 1/c correction,
larger for higher zeros (at c = 28.2 the deficit is 0.21, 0.42, 0.54,
0.68 for zeros 1 to 4), so higher zeros need larger windows to settle.
The sign of the deficit shows it comes from the finite-window vacuum
differing from Phi, not from the shares w_gamma (which could only raise
K for the first zero).

**Where this leaves the response theory.** Everything in it is now
classical except one function of the window:

    tick_n(c) = gamma_n + lambda_0(c) / (gamma_n Xi~'(gamma_n)^2) x (1 + O(1/c))

and a mirror at lag a changes lambda_0 by the self-overlap of Phi, with
the level/shape cancellation rho_n(a) also from Xi alone. The only thing
the construction knows that the classical function does not is
lambda_0(c): the vacuum's level, and the fact that it is positive.

## The one curve: the vacuum's level against the window (`level_curve.py`)

lambda_0(c) for both halves, c = 4.00 to 9.00 in steps of 0.05, fixed
N = 80 (the earlier stored runs mixed basis sizes, which shows up as a
2 to 5% zigzag; at fixed N the truncation error is 0.005 to 0.03 decade
and smooth in c). Data: `level_curve.json`.

The level falls from 1e-12 to 1e-37 (even) over this range. Its rate of
fall, in decades per unit of c (4 pi / ln 10 = 5.457 is the asymptotic
rate):

- mean 5.08 (even), 4.93 (odd), so still below the asymptotic rate;
- **it is not steady**: it swings between 4.4 and 5.8, rms 0.37, with
  bumps roughly 0.4 to 0.6 apart, the two halves partly out of step
  (the mean of the two has rms 0.28);
- **kinks where a prime power enters the window**: the rate jumps by
  +0.41 / +0.68 (even / odd) at c = 5, +0.84 / +1.07 at c = 7,
  +0.26 / +0.18 at c = 8, and nothing at c = 6 (-0.17 / -0.07). A new
  mirror makes the level fall faster from the moment it enters.

**Tested and not found:** a tie between the bumps and the zeros. Binning
the rate by where 2 pi c sits between neighbouring zeros gives 4.93,
5.03, 5.19, 5.15 (even, +-0.08) and 5.01, 4.81, 4.93, 5.00 (odd, +-0.07);
the controls with pi c and 4 pi c look the same. The mean spacing of
zeros in c at this height (about 0.5) matches the spacing of the bumps,
but the bumps do not line up with the zeros.

**Not understood:** what the bumps are, and the size of the kinks. A
first-order estimate (a mirror entering at the edge adds -2 w beta^2 (L -
ln p) to the level, beta the vacuum's edge value) would multiply the rate
by 1 + ln p / sqrt p = 1.7 at c = 5 and 7; the measured jump is 8 to 25%.
Either the edge bookkeeping (Hadamard's formula for this form) is not
what was assumed, or the vacuum rearranges at once. Next: measure the
edge value directly and test d(lambda_0)/dL against it.

### A prime entering the window (`prime_entry.py`, prime 7, N = 120)

E = w <v0|A_ln p|v0> / lambda_0: the mirror's term in the vacuum's level,
in units of the level. By Hellmann-Feynman, E = - d ln(lambda_0) / d ln(w)
exactly. h = ln(c / 7) is the width of the strip the mirror connects.

| c | 7.005 | 7.01 | 7.02 | 7.04 | 7.08 | 7.15 | 7.3 |
|---|---|---|---|---|---|---|---|
| h | 0.0007 | 0.0014 | 0.0029 | 0.0057 | 0.0114 | 0.0212 | 0.0420 |
| E | 0.030 | 0.078 | 0.28 | 2.0 | 35 | 1240 | 262,000 |

- **Entry is smooth.** Just inside (0.07% past the prime) the mirror
  carries 3% of the level. There is no jump, which is why the measured
  kink is small.
- **Then its share explodes.** Neither the naive estimate (E growing like
  h) nor the log-softened one fits: E / h rises from 42 to 6e6. This is
  the vacuum's steep rise away from the edge (the Bessel profile
  I_0(C sqrt t) of `string.md`, C about sqrt 2 x 2 pi c): the strip the
  mirror connects soon reaches where the vacuum is exponentially larger
  than at the edge.
- **So the level is a remainder.** At c = 7.3 the level is
  (everything else) - (the mirror at 7), two terms each 262,000 times
  larger than their difference. One relative part in a million of the
  prime's weight would change the level by a quarter. For older primes
  the same ratio is w x 2 R(ln p) / lambda_0, the "strength in units of
  the vacuum's depth" of the fade-in section: tens of decades.
- E reaches order 1 at h of about 0.004, i.e. c about 7 (1 + 0.004):
  the same place where the form without the prime would go negative
  (the margin 1 / (8 p^2) = 0.0026 of `string.md`).

So the two readings offered before were both wrong. The rate of fall
changes little at entry because the mirror starts with nothing to
connect; within a few percent of window growth it has become one side of
a cancellation to five or six digits, and the vacuum keeps adjusting so
that the remainder stays positive and keeps falling at about the same
rate.

### The bumps in the level's rate of fall: three tests (`bumps.py`, `level_curve.json`)

1. **New near-null states arriving (once per unit of c in each half):
   no.** Folding the detrended rate by the fractional part of c explains
   4% (even) and 5% (odd) of its variance. The bumps do not repeat with
   period 1.
2. **Are they an artifact of the basis size: no.** Recomputed with
   N = 120 for c = 5.5 to 6.55: the level shifts by a slowly varying
   0.01 to 0.03 decade, and the rates agree with N = 80 to +-0.07 while
   swinging between 4.5 and 5.7. The bumps are real.
3. **Zeros, looked at again: suggestive, not established.** Averaging the
   detrended rate at c = gamma / 2 pi + s over the 8 zeros in range:
   even half peaks at +0.17 for s = -0.15 (smooth profile: -0.15, -0.04,
   +0.10, +0.17, +0.14, +0.07, +0.06, -0.06 for s = -0.30 ... +0.05); odd
   half peaks at +0.13 for s = 0. So the even rate tends to be highest
   shortly before 2 pi c reaches a zero and the odd rate when it gets
   there, which is also why the earlier four-bin test saw little. But a
   randomly displaced comb does as well in 12% (even) and 22% (odd) of
   trials, and the peaks are half the size of the bumps (rms 0.35).
   Eight zeros are too few. The maxima not near a zero sit just after
   the prime 7 and near c = 4.5 and 8.6.

The periodogram's "period 0.67" is not a clean period: successive maxima
are 0.3 to 0.6 apart, irregular, like the zeros.

To decide: extend the curve to c of about 16 (some 25 more zeros).

### The curve extended to c = 16, and the zeros test decided (`level_curve.py` N = 120, 141 windows from 9 to 16)

The level falls from 1e-37.5 to 1e-74.4 (even half). Rate of fall: mean
5.27 (even), 5.19 (odd), closer to 5.46 than at smaller windows; range
4.5 to 6.1; rms about its trend 0.26 / 0.25 (it was 0.37 for c = 4 to 9).

**Zeros: no.** Out-of-sample test with the shifts fixed from the first
data set (even: s = -0.15, odd: s = 0), 16 new zeros: zero-triggered mean
of the detrended rate -0.001 (z = 0.0) and -0.020 (z = -0.5). The whole
profile for s between -0.3 and +0.3 stays within +-0.09 with no shape.
The pattern seen with 8 zeros was chance. The bumps are not tied to the
zeros crossing 2 pi c.

**Primes: yes, cleanly.** Jump of the rate across an integer (after minus
before, even / odd):

| c | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|
| jump | +0.04 / -0.22 | +1.11 / +0.91 | 0.00 / +0.07 | +1.33 / +0.96 | -0.22 / +0.21 | +0.04 / +0.14 |

Together with 5 (+0.41 / +0.68), 7 (+0.84 / +1.07), 8 (+0.26 / +0.18)
and 6 (none): every prime gives a jump of about one decade per unit, a
prime power a smaller one, other integers none.

**What the bumps look like now.** After a prime the rate is high (5.8 to
6.1) and relaxes over roughly a unit of c to about 5.0, reaching its
lowest values just before the next prime (10.97: 4.71 / 4.45; 12.97:
4.77 / 4.47; 6.97: 4.87 / 4.28; 4.97: 4.88 / 4.25). So part of the
"bumps" is a sawtooth driven by the primes: a kick at each entry, then a
decay. That also explains why the swings were larger between 4 and 9,
where prime powers arrive more often. On top of the sawtooth there
remain wiggles of about +-0.3, 0.4 to 0.6 apart, in both halves, weakly
correlated (+0.21), not yet attributed to anything.

Ruled out so far for those wiggles: the arrival of near-null states
(period 1), the basis size, the zeros at 2 pi c. Not yet tested: the
vacuum's edge value itself, and the powers of 2 and 3 inside the window.

### The wiggles: a wave through the ladder (stored spectra `ladder_fine.json`, c = 5.00 to 7.62, step 0.02, N = 80)

Looking at the levels above the vacuum, not only at the vacuum:

- **Neighbouring levels wiggle against each other.** Correlation of the
  rates of fall: level 0 with level 1: -0.38 (both halves); 1 with 2:
  -0.26 / -0.23; 0 with 2: +0.36 / +0.44. The sign alternates up the
  ladder.
- **The gaps between successive levels all turn at the same windows.**
  Zero crossings of the detrended gaps, even half: levels 0-1 at 5.03,
  5.39, 5.79, 6.16, 6.42, 6.65, 7.10, 7.41; levels 1-2 at 5.12, 5.38,
  5.68, 6.10, 6.45, 6.68, 7.04, 7.39; levels 2-3 at 5.13, 5.50, 5.68,
  5.99, 6.44, 6.73, 7.05, 7.39. Half-period about 0.33, so a full period
  of about two thirds of a unit of c (the 0.67 of the periodogram, and
  the spacing still seen between c = 13 and 16).
- **The two halves are a quarter period out of step.** The odd gaps cross
  zero 0.15 to 0.2 later than the even ones; correlation of the lowest
  even and odd gaps: +0.21 unshifted, -0.64 and +0.55 with a shift of
  +-0.25 in c.
- **Averaged over the ladder the wiggle fades.** Rms of the detrended
  mean log of the lowest k levels: 0.032, 0.014, 0.011, 0.008 for
  k = 1, 2, 3, 4 (even half).

So the wiggle is one wave running through the whole ladder: each step up
the combined ladder (even, odd, even, ...) shifts its phase by a quarter
turn. It is consistent with every level carrying the same wiggle as a
function of its "age" (window size minus the size at which it was born),
since levels are born half a unit apart and 3 pi x 0.5 is a quarter turn
backwards. The levels trade a little with their neighbours while the
family as a whole (the product of the tiny levels) falls smoothly.

Not explained: why the period is about 2/3 of a unit of c, and whether
it is exactly that.

Summary of the fine structure of the vacuum's level: a smooth fall
(approaching 4 pi per unit of c) + a kick of about one decade per unit at
each prime, decaying over about a unit + a wiggle of about +-0.3 that is
the vacuum's share of a wave through the ladder. Nothing from the zeros.

### "Why 2/3?": the wiggle is a set of sharp tones (periodograms of `level_curve.json`)

Rate of fall of the vacuum's level, quadratic trend removed, windows
within 0.12 of a prime power left out. Strongest frequencies in cycles
per unit of c (share of variance in brackets):

| data | even half | odd half |
|---|---|---|
| c = 9 to 16, N = 120 | 1.495 (1.0), 2.630 (0.33), 1.705, 3.760, 1.285 | 1.495 (0.87), 2.625 (0.41), 0.845, 1.710, 1.275 |
| c = 4 to 9, N = 80 | 1.485 (1.0), 2.640 (0.45), 3.740, 0.790, 2.165 | 1.505 (1.0), 2.655 (0.40), 0.780, 3.735, 2.110 |

- **The main tone is at 1.50 cycles per unit, in both halves and in both
  independent ranges** (1.485 to 1.505; the resolution is 0.14 to 0.2 but
  a strong line is located much better than that). So the period is 2/3
  to within about half a percent, and the phase is 3 pi c = 2 pi c x 3/2.
- **There are overtones**, the same in both ranges: 2.64 +- 0.01, 3.74,
  and weaker ones near 0.79, 2.13, 1.28, 1.71. They are not multiples of
  1.5.
- Candidate readings, none established: simple fractions (3/2, 8/3 =
  2.667, 15/4, 4/5, 9/7, 12/7), which fit 2.64 only at two standard
  deviations; or pairs differing by sqrt 2 (1.5 and 2.12; 2.646 = sqrt 7
  and 3.742 = sqrt 14), which fit the numbers better but have no
  explanation.

A quarter turn per ladder step is the Fourier eigenvalue i^n of these
waves, so the wiggle looks like an interference between each wave and
its own tail outside the window (the tail of a function cut off at
lambda oscillates as sin(2 pi lambda x - n pi / 2)); a phase 3 pi c would
be that tail read at x = 3 lambda / 2, i.e. the ratio 3/2 of the first
two primes. That is a guess. To separate fractions from square roots the
curve has to be extended (to c = 30 the resolution would be 0.04).

### The curve extended to c = 24, and the tones again (`level_curve.py` N = 150, 80 windows from 16.1 to 24.0, step 0.1)

The level reaches 1e-117 (even) and 1e-112.5 (odd) at c = 24. Mean rate
of fall 5.32 / 5.28. Prime kicks again: +0.83 / +0.99 at 17, +1.15 / +1.01
at 19, +1.33 / +1.59 at 23; at 18, 20, 21, 22 between -0.5 and +0.3.

**Cleanest view of the wiggle: the difference of the two halves' levels**,
log10(even) - log10(odd), cubic trend removed. The prime kicks are nearly
the same in both halves and cancel. Strongest lines by range:

| range of c | rms of the wiggle | lines (cycles per unit of c; share) |
|---|---|---|
| 4 to 9 | 0.041 | 1.500 (1.1), 0.772, 2.630 |
| 9 to 12.5 | 0.030 | 1.496 (1.2), 2.606, 2.160 |
| 12.5 to 16 | 0.026 | 1.470 (1.4), 1.038 |
| 16.1 to 20 | 0.023 | 1.486 (1.1), 0.894 (0.5) |
| 20 to 24 | 0.016 | 0.898 (0.8), 1.248; no line at 1.5 |

- **The main tone sits at 1.50 in four successive ranges** (1.500, 1.496,
  1.470, 1.486; each range resolves about 0.25, a strong line is placed
  to a few hundredths). It does not drift. So the pitch is 3/2 as far as
  this can tell.
- **It fades.** The wiggle's size falls from 0.041 to 0.016, and beyond
  c = 20 the 3/2 tone is gone; a slower tone near 0.9 (0.894, 0.898)
  dominates there. The overtone near 2.63 is present up to c of about 12
  and not after.
- In the rates of the single halves (previous section's method) the
  kicks at 17, 19, 23 leak into the spectrum and made the lines look as
  if they drifted (1.44 to 1.47, 2.39); the difference of the halves
  does not show that.

So "why 2/3" becomes: why a tone at exactly 3/2 that dies out as the
window grows, and what replaces it. A tone that fades fits something
tied to particular small primes (their influence relative to the whole
shrinks as more primes arrive) better than something built into the
waves themselves. The guess "ratio of the first two primes" survives;
the square-root reading does not need to be kept, since the 2.63 line
was never accurate enough and has gone. Next candidates to check against
the later tone near 0.9: 8/9 = 0.889 and 9/10; needs c beyond 24.

### Testing the reading of the 3/2 tone (fits to `level_curve.json`)

The "2 and 3" idea in its concrete form: each ladder wave interferes
with its own tail outside the window, read at x = 3 lambda / 2. A
function cut off at lambda has a tail sin(2 pi lambda x - n pi / 2) / x,
so this predicts three things at once: phase 2 pi lambda (3 lambda / 2) =
3 pi c; amplitude falling like 1 / x, i.e. like 1 / lambda = c^(-1/2);
and a quarter turn between successive ladder levels (even and odd
vacua).

Fit of the detrended difference of the halves' levels to
c^(-p) (a cos 2 pi f c + b sin 2 pi f c), everything free:

| range | frequency | amplitude | variance explained |
|---|---|---|---|
| 4 to 20 | 1.4969 | 0.113 c^(-0.52) | 0.59 |
| 4 to 12 | 1.5022 | 0.129 c^(-0.59) | 0.58 |
| 12 to 20 | 1.4929 | (flat within this short range) | 0.62 |

With the frequency fixed at 3/2 and the amplitude at c^(-1/2), each
half's own level (c = 9 to 20): even 0.088 / sqrt(c) at phase -165
degrees, odd 0.084 / sqrt(c) at +118 degrees. The difference is 283
degrees, i.e. the odd one is 77 degrees behind: a quarter turn within
the accuracy.

**All three signatures hold**: frequency 1.497 (3/2 to 0.2%), amplitude
exponent -0.52 (predicted -0.5), relative phase 77 degrees (predicted
90). A single such tone accounts for about 60% of the wiggle.

What this does and does not show. It supports "the tail at one and a
half window lengths". It does not tell whether 3/2 is there as the ratio
of the primes 3 and 2 or as the integers 3 and 2 (in the sum over
integer multiples that builds these states, three times the point
lambda / 2 where twice the point reaches the edge). Tones at other simple
ratios (4/3, 5/3, 2, 5/2, 3) are not seen, so neither a "pairs of
primes" nor a "pairs of integers" family is visible; 3/2 stands alone.
The decisive experiment is a form for a number system without the prime
2 or without the prime 3 (a Dirichlet L-function mod 4 or mod 3): if the
tone is about the primes it must vanish or move there; if it is about
the integers it stays. That needs the form for L-functions, which this
code does not have yet.

### Ludvig's point on the missing cousins of 3/2 (27 Nov, 05:20)

"Are the tones not just ratios of the prime pairs?" This removes the
objection that 5/3, 5/2, 7/2 and so on are not seen. If a tone at p/q is
weighted by how strongly both primes act on the vacuum, the fade-in
kernel says only 2 and 3 sit in the bright part: 5 acts with about 1e-6
of the strength of 2, 7 with 1e-12. So 3/2 should be the only audible
ratio, and every pair involving 5 or more should be at least a million
times fainter. The absence of the cousins is what the prime-pair reading
predicts, not evidence against it. (The same weighting makes a
"pairs of integers" family nearly silent too, apart from 3/2 and the
powers of 2 and 3, so this does not by itself separate the two.) The
faint lines near 2.63, 0.8 to 0.9 and 3.74 are not ratios of small
primes and are left unexplained. Decision: park the tone; it will be
checked for free when the form for an L-function (no prime 2, or no
prime 3) is built for other reasons.

## Another number system: the form for a Dirichlet L-function (`lfunc_block.py`, `lticks.py`)

Built from the explicit formula for a real primitive character chi of
conductor Nc and parity a, as a change to zeta's kernel function H:

- each prime power q carries chi(q); primes dividing Nc drop out;
- no pole: zeta's term -xU (measure -dt) is replaced through the trivial
  zeros by -U (a = 0) or -xU - U + (x + 1) log((x + 1)/2) (a = 1); the
  difference of measures is (1 + 1/t) dt, resp. (1 + 1/(t(t + 1))) dt,
  which is smooth, so the regularisation at t = 1 is unchanged;
- the conductor adds log(Nc) g(0), i.e. -(1/2) log(Nc)(x - 1) in the
  bracket.

**It worked at the first run.** c = 24.2, N = 90:

| | chi mod 4 (no prime 2) | chi mod 3 (no prime 3) |
|---|---|---|
| lowest even eigenvalues | 6.8e-29, 1.9e-21, 1.9e-15, 3.4e-10 | 1.7e-39, 8.7e-32, 6.0e-25, 2.2e-19 |
| lowest odd eigenvalues | 5.4e-25, 2.9e-18, 8.1e-13, 4.8e-8 | 1.7e-35, 2.8e-28, 5.2e-22, 1.0e-16 |
| ticks of the lowest even state | 6.020948905, 10.2437703, 12.98809801, 16.3426071, 18.2919932, 21.45061134 | 8.039737156, 11.24920621, 15.70461918, 18.2619975, 20.45577081 |
| tick minus zero of L(1/2 + it, chi) | 3e-27, 3e-25, 8e-24, 6e-22, 6e-21, 2e-19 | 9e-37, 3e-35, 7e-33, 3e-31, 3e-30 |

So the form is positive, has the ladder of tiny levels, and the recipe
gives the zeros of the L-function (computed independently with mpmath)
to 27 and 37 digits. Tick minus zero is again a multiple of the lowest
eigenvalue (47 and 500 for the first zero).

Fewer tiny levels than for zeta at the same window: the zeros are denser
by log(Nc), so the "sea" reaches only 2 pi c / Nc, and there are about
2c / Nc near-null states (about 12 for chi mod 4 at c = 24.2, 16 for chi
mod 3).

### The level curves of the two other systems, and what the tones are (`CHI=chi4|chi3 level_curve.py`, 201 windows each, c = 10 to 30, N = 110; `level_curve_chi.json`)

Mean rate of fall: 1.32 / 1.27 decades per unit of c for chi mod 4 and
1.77 / 1.72 for chi mod 3, against 4 pi / (Nc ln 10) = 1.364 and 1.819.
So the natural window is c / Nc, as the count of tiny levels said.

**Kicks follow the character.** Jump of the rate at a prime (even / odd):
chi mod 4: +0.04 / +0.24 at 13, +0.12 / +0.20 at 17, +0.14 / +0.26 at 29
(chi = +1); -0.03 / -0.25 at 11, -0.16 / -0.20 at 19, -0.18 / -0.27 at 23
(chi = -1). chi mod 3: +0.20 / +0.21 at 13, +0.17 / +0.33 at 19 (chi = +1);
+0.03 / -0.36 at 11, -0.17 / -0.26 at 17, -0.32 / -0.32 at 23,
-0.30 / -0.38 at 29 (chi = -1). A prime with chi = +1 speeds the fall, one
with chi = -1 slows it. In units of c / Nc the size is about one decade
per unit, as for zeta. Non-prime-powers: within +-0.1.

**Tones** of the detrended difference of the halves' levels, converted to
cycles per unit of c / Nc:

| system | primes it has | lines |
|---|---|---|
| zeta | all | 1.50, 2.63 to 2.65, 3.74, 0.77 to 0.90, 2.13, 1.04, 1.71, 1.28 |
| chi mod 3 | 2, 5, 7, ... (no 3) | **1.500** (share 1.35), 2.12, 1.71, 3.74, 1.30, 0.52. Nothing at 2.65. |
| chi mod 4 | 3, 5, 7, ... (no 2) | **2.648** (share 1.02), 1.06, 4.77, 2.91, 0.64, 0.35. Nothing at 1.5 (share 0.001). |

**So the 1.5 tone needs the prime 2 and not the prime 3; the 2.65 tone
needs the prime 3 and not the prime 2.** "3/2" was never the ratio of
two primes.

**The rule that fits: a mirror at q gives a tone at q - 1/q.**

| q | q - 1/q | seen at | in which systems |
|---|---|---|---|
| 2 | 1.500 | 1.497, 1.500 | zeta, chi3 (have 2); absent in chi4 |
| 3 | 2.667 | 2.63 to 2.65 | zeta, chi4 (have 3); absent in chi3 |
| 4 | 3.750 | 3.74 | zeta, chi3 (have 4); absent in chi4 |
| 5 | 4.800 | 4.77 | chi4 (zeta not scanned that high) |
| 3/2 | 0.833 | 0.77 to 0.90 | zeta only (needs 2 and 3) |
| 5/2 | 2.100 | 2.12, 2.13 | zeta, chi3 (have 2 and 5); absent in chi4 |
| 5/3 | 1.067 | 1.04, 1.06 | zeta, chi4 (have 3 and 5); absent in chi3 |

Eleven line-and-system matches and seven predicted absences, all correct.
The measured values sit 0.2 to 0.7% below q - 1/q. Unidentified lines:
1.71 and 1.28 to 1.30 (zeta, chi3), 0.52 (chi3), 2.91, 0.64, 0.35 (chi4).

Reading: q - 1/q is the difference of q and 1/q, the two images of the
window's edge under the mirror (x = q lambda and x = lambda / q). The
phase 2 pi lambda (q lambda - lambda / q) is the beat between the wave's
tail read at those two points. So each prime in the bright part of the
vacuum plays its own note, and pairs of them play the notes of their
ratios. Ludvig's question ("are the tones not just ratios of the prime
pairs?") was half of it: the pair tones are there (3/2, 5/2, 5/3), and
the loudest line is the single prime 2.

### Passivity (Ludvig, 27 Nov 06:00): an ill-designed test, and where the idea already lives (`interlace.py`)

Suggestion: positivity as energy conservation, the system as passive.
Test tried: do the zeros of the even vacuum's cosine transform and of the
odd vacuum's sine transform interlace, as resonances and anti-resonances
of a lossless passive system do? Result at c = 8.2: the two sets are the
same thirteen numbers, the zeros of zeta, equal to 5.6e-28. Of course:
both halves of the form are sums over the same zeros, so both vacua
vanish at all of them. The test carried no information; it was the wrong
pair of functions. (The run at c = 14.2 with default quadrature is below
its precision floor and is discarded.)

Where passivity already is, in the literature this rests on: for every
finite window the recipe's operator is self-adjoint, so its ticks are
real whatever the zeros of zeta do (Connes-van Suijlekom; the hypothesis
is that the lowest eigenvalue is simple and even). Each window is a
lossless system by construction. The Riemann hypothesis would follow if
the ticks converge to the zeros of zeta as the window grows (Hurwitz),
and that convergence is the open step. So "find a passive system" is
done window by window; what is missing is that the sequence of passive
systems has zeta as its limit. Burnol's "causality is equivalent to RH"
and de Branges' strings are the same idea in other clothes
(`string.md`).

## A bound to aim for: the distance of a tick from its limit, from the window's own data (`tick_bound.py`)

For the lowest even state v (unit norm, coefficients v_m on the
frequencies w_m), its cosine transform v^, its n-th tick t_n and the
lowest eigenvalue eps, define, with nothing but the form on the window c,

    B_n(c) = t_n P_n eps / v^'(t_n)^2 ,     P_n = sum_m v_m^2 / (t_n^2 - w_m^2).

Pairing the eigen-equation with v divided by (s^2 - t_n^2) gives
t_n - gamma_n = B_n (1 + r_n), where r_n = 2 sum over the other zeros of
w_gamma / ((gamma^2 - t_n^2) P_n) and the w_gamma are the shares of eps
carried by each zero. (The term with cos(t_n u) drops out of P_n because
v^(t_n) = 0.)

**Measured rho_n = (t_n - gamma_n) / B_n**, true zeros from mpmath:

| system, c | eps | tick 1 | tick 2 | tick 3 | tick 4 |
|---|---|---|---|---|---|
| zeta, 4.2 | 1e-13.0 | 1.0501 | 1.1285 | 1.2003 | 1.3521 |
| zeta, 6.2 | 1e-23.0 | 1.0204 | 1.0488 | 1.0718 | 1.1126 |
| zeta, 8.2 | 1e-33.3 | 1.0100 | 1.0236 | 1.0342 | 1.0523 |
| zeta, 10.2 | 1e-43.8 | 1.0066 | 1.0154 | 1.0222 | 1.0336 |
| zeta, 12.2 | 1e-54.3 | 1.0045 | 1.0104 | 1.0149 | 1.0224 |
| zeta, 14.2 | 1e-64.8 | 1.0032 | 1.0074 | 1.0106 | 1.0159 |
| zeta, 16.2 | 1e-75.4 | 1.0024 | 1.0056 | 1.0080 | 1.0119 |
| chi mod 4, 14.2 | 1e-15.0 | 1.0095 | 1.0306 | 1.0513 | 1.0858 |
| chi mod 4, 26.2 | 1e-30.8 | 1.0027 | 1.0085 | 1.0139 | 1.0224 |
| chi mod 3, 14.2 | 1e-21.1 | 1.0087 | 1.0181 | 1.0370 | 1.0514 |
| chi mod 3, 26.2 | 1e-42.3 | 1.0026 | 1.0053 | 1.0106 | 1.0145 |

(also zeta at 5.2, 7.2, 9.2 and both characters at 20.2: same pattern.)

- **In every case rho_n is above 1 and tends to 1.** 64 values, three
  systems.
- **rho_n - 1 = 0.0032 x (gamma_n Nc / c)^2 to within 10%**, for every
  tick, window and system: (rho_1 - 1) c^2 = 0.88, 0.71, 0.78, 0.64,
  0.67, 0.65, 0.69, 0.67, 0.65, 0.63 for zeta at c = 4.2 ... 16.2; the
  higher ticks scale with gamma_n^2 (tick 2 / tick 1: 2.33 against 2.21;
  tick 4 / tick 1: 4.96 against 4.63); the characters give 0.0032 and
  0.0030 with c replaced by c / Nc.
- Reading: r_n is the average of gamma_n^2 / gamma^2 over the zeros that
  carry the vacuum's level, so those zeros sit around
  gamma = 17.7 c / Nc = 2.8 x (2 pi c / Nc): just beyond the reach of the
  window, as expected.

**The conjecture this suggests** (a target, not a result):

    B_n(c)  <=  t_n(c) - gamma_n  <=  B_n(c) (1 + C (gamma_n Nc / c)^2),   C about 0.0035,

for every window. Its content: the distance from a tick to its limit is
given, to a relative accuracy that improves like 1/c^2, by a quantity
computed on the window alone. If something of this shape could be proved
from the form's side, without assuming where the zeros are, it would be a
rate of convergence for the recipe. As it stands the identity behind it
assumes the form is a sum over real zeros.

Practical by-product: t_n - B_n (1 + r_n) improves each tick by two to
three further digits.
