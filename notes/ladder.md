# The top of the ladder: states one by one as the window grows

Started 26 Nov 2026, late. Aim : explore the few
shallow states where the known theory runs out (Andrews-Connes-La Gatta
reduce positivity to "transition modes at the edge of the ladder"), with
their open problem "explain mu_lambda about 2/lambda" as the target.

## Step 1: the terrain

`scripts/ladder_track.py 100 10 c...`, data
`ladder_track.json`: the ten lowest eigenvalues of each half of the form,
c = 1.5 to 8.5 in steps of 0.1, N = 100, ball arithmetic.

**Births happen at integers and half-integers, not at primes.** A level
"is born" when it drops below 0.1 on its way down.

| half | window c at which successive levels pass 0.1 |
|---|---|
| odd | 1.96, 2.97, 3.98, 4.98, 6.00, 6.99, 8.01 |
| even | 2.45, 3.48, 4.46, 5.49, 6.45, 7.50, 8.48 |

- One new tiny level per unit of c in each half, the even half half a
  unit ahead. 4, 6 and 8 get a birth exactly like 2, 3, 5 and 7. So the
  plunge is phase-space counting (the gas in the box: one more mode fits),
  not arithmetic. This agrees with Connes-Consani's remark that the m-th
  eigenvector matches once mu > m (odd: mu > m + 1) and with their
  nu(mu) = 2 mu - 1 at half-integers.
- Once plunging, a level falls by about one decade per 0.2 in c, steadily
  (4 pi per unit c in the exponent, less the power-law prefactor).

**Before it plunges, a level comes down a staircase.** The seventh odd
level sits near 0.89 (c = 5.4 to 6.2), then near 0.74 (c = 6.6 to 7.0),
then starts its plunge at about 7.2. Other levels show the same plateau
values (0.89, about 0.85, 0.74 to 0.78) at other windows. The steps
between plateaus are where to look for the primes; step 0.1 in c is too
coarse to say where they are.

**Bottom of the ordinary spectrum** (first eigenvalue above 0.3) at
integer c, times sqrt(c)/2, which would be 1 if it were 2/lambda:

| c | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---|---|---|---|---|---|---|
| odd | 0.66 | 0.72 | 0.84 | 0.87 | 0.91 | 0.96 | 0.94 |
| even | 0.47 | 0.60 | 0.55 | 0.68 | 0.64 | 0.70 | 0.65 |

The odd half approaches the 2/lambda rule from below. The even half is
lower and jumpy, as expected here: Andrews-Connes-La Gatta's number is
for the form with the pole direction also removed, which we have not
done. So this table is not yet their quantity.

**Reading.** The tiny levels are the gas; they do not care which integers
are prime. What the primes do must be visible elsewhere: in the margin
(margin.md), and possibly in the staircase of the level that is next in
line.

## Step 2: the staircase, zoomed (odd half, c = 5 to 7.6)

Data: `ladder_fine.json` (step 0.02, N = 80), `ladder_zoom.json` (step
0.005 around 6 and 7), `ladder_no7.json` (the same with the term of 7
left out; `EXCL=7`).

**Three kinds of state at the top of the ladder.**

1. **Flat states.** Levels that do not move as the window grows: 0.8909
   (c = 5.60 to 5.72), 0.7741 (5.64 to 5.72), 0.8771 (6.12 to 6.24),
   0.7369 (7.12 to 7.24), each constant to four digits over those ranges.
   By the Hadamard formula d eps/da = -2 beta^2, flat means beta = 0: the
   state has no amplitude at the window's ends. They are also untouched
   by the prime: with or without the term of 7 the flat level is 0.7369
   against 0.7366 at c = 7.1 to 7.3.
2. **The faller.** One state at a time comes down through the flat ones
   by narrow avoided crossings (near 7.00: gap 0.015, the faller's slope
   -0.74 per unit c before, the flat one's 0), which is what made the
   "staircase" of step 1. After its last crossing it becomes the next
   tiny level. There is no sharp kink at a prime; at c = 6 (no term)
   everything is smooth as well.
3. **Without the prime: a diver.** Leave out the term of 7 and, just past
   7, one eigenvalue goes negative and then falls explosively:
   -6e-25 (7.03), -2e-18 (7.04), -2e-13 (7.05), -7e-10 (7.06), -1e-7
   (7.08), -8e-5 (7.10), -0.0085 (7.15), -0.042 (7.20), -0.107 (7.30).
   (N = 80 does not resolve the margin itself, which is 7.0175.) The
   faller also comes down faster without the prime: 0.655 against 0.679
   at c = 7.1.

**What the prime does, seen from here.** It does nothing to the flat
states, brakes the faller mildly, and prevents the diver. The diver starts
in the deepest states (that is the margin of margin.md) but within 0.3 of
window it is an order-one negative level: the damage does not stay small.

**The 2/lambda rule, first reading.** The bottom of the ordinary spectrum
is, at any window, the lowest flat level not yet passed by the faller:
0.850 near c = 4.9, 0.774 near 5.7, 0.737 near 7.2, against 2/sqrt(c) =
0.90, 0.84, 0.745. So the smooth rule looks like the envelope of a
descending staircase of flat (edge-blind) levels. If so, "explain
2/lambda" becomes "what are the flat states and why are their levels
spaced as they are". Not yet examined: their shapes, the even half, and
the form with the pole direction removed (which is the quantity in
Andrews-Connes-La Gatta).

**Physics words for the flat states (ours, tentative).** A state that has
no amplitude where the system couples to the outside does not decay and
does not shift: a dark state, or a bound state in the continuum. Here the
"outside" is reached through the window's ends, and both the window's
growth and the primes' mirrors act only there.

## Step 3: what the flat states are (odd half, c = 7.2, N = 360)

`deep_profile.py 360 7.2 9` (now also saves the eigenvectors,
`deep_profile_7.2_360.npz`); analysis inline, numbers below.

**Check: each eigenvalue is the state's weight on zeta's zeros.** With
fhat the transform in scaling frequency, 2 sum over zeros of
|fhat(gamma)|^2 (first 59 zeros) gives 0.6295, 0.7370, 0.889 for levels
7, 8, 9 against eigenvalues 0.6430, 0.7371, 0.8935. So "low eigenvalue"
means "little of the state's spectrum lands on a zero".

**Where the states sit in the window** (share of weight by distance from
the nearest end; half-window 0.987):

| level | eps | 0 to 0.02 | 0.02 to 0.1 | 0.1 to 0.4 | 0.4 to 0.99 |
|---|---|---|---|---|---|
| 1 | 1.5e-25 | 0 | 0 | 0.0004 | 0.9996 |
| 3 | 4.1e-14 | 0 | 0 | 0.715 | 0.285 |
| 5 | 2.2e-5 | 0.0003 | 0.180 | 0.641 | 0.179 |
| 6 (newest tiny) | 0.026 | 0.047 | 0.536 | 0.135 | 0.283 |
| 7 (faller) | 0.643 | 0.055 | 0.094 | 0.281 | 0.570 |
| 8 (flat) | 0.737 | 0.0006 | 0.034 | 0.591 | 0.374 |
| 9 (flat) | 0.894 | 0.0048 | 0.005 | 0.254 | 0.736 |

- The tiny states are nested shells: the deepest in the middle, each
  later one nearer the ends, the newest hugging them.
- The flat states keep away from the ends (0.06% and 0.5% of their
  weight in the outer 2%), the faller does not (5.5%). This is the
  "no amplitude at the ends" inferred in step 2, seen directly.

**What the flat state looks like in frequency.** Level 8 is a comb
between frequencies 44 and 64: peaks at 46.3, 50.8, 54.8, 58.5, 62.8, and
its transform vanishes at 49.0, 52.75, 56.75, 60.0, 66.0. Zeta's zeros
there: 48.0, 49.8, 53.0, 56.4, 59.3, 60.8, 65.1, 67.1. The comb's gaps
sit on single zeros (52.75 on 53.0, 56.75 on 56.4) and between the two
members of a close pair (49.0 between 48.0 and 49.8; 60.0 between 59.3
and 60.8; 66.0 between 65.1 and 67.1). The state dodges the zeros as well
as a window of this length allows; 0.737 is what it cannot dodge.

**Counting (known, Connes 1999 in other words).** A window of length
L = ln c has c ln c sine modes and c ln c cosine modes below the
frequency W = 2 pi c. Zeta has 2(c ln c - c) zeros of size below W. The
difference is 2c: the number of tiny states. Below W there are more
modes than zeros, so a state can vanish at every zero: the tiny states.
W is where the zeros' spacing 2 pi / ln(omega / 2 pi) equals the
window's resolution 2 pi / L. Above W the zeros outnumber the modes and
not all can be dodged: the ordinary spectrum.

**A candidate reason for the 2/lambda rule (heuristic, ours, untested).**
Just above W the excess of zeros over modes in a band of width D is about
D^2 / (4 pi W). It reaches one zero at D = sqrt(4 pi W) = sqrt(8 pi^2 c),
which is 23.8 at c = 7.2; the flat state's band is 44 to 64, width 20.
So the lowest ordinary state spreads over the widest band in which it
still has (almost) as many modes as zeros, and that band's width grows
like sqrt(c) = lambda. The weight it must leave on the one zero it cannot
dodge is of order one over the number of modes in the band, D L / 2 pi,
i.e. of order 1 / (lambda ln c) up to constants. That gives the 1/lambda
of mu = 2/lambda, with a logarithm the measured rule does not show; to be
tested against the level's dependence on c before believing it.

**In the vacuum words of picture.md:** tiny states = modes that
fit entirely between the absorption lines; flat states = the lowest
excitations above the Fermi frequency, kept away from the walls; their
level = a zero-point cost set by one unavoidable zero.

## Step 4: the heuristic tested on c = 5 to 30. It fails; what replaces it is sharper

`flat_track.py 300 c...` (odd half; the form in ball arithmetic, then
double-precision eigenvectors, which is enough for eigenvalues above
0.2), data `flat_track.json`; zeros in `zeta_zeros_139.json`.

**The prediction of step 3 is wrong.** The low ordinary states do not sit
in a band of width sqrt(8 pi^2 c) just above W = 2 pi c. From c = 13 on,
the lowest ones are narrow in frequency (a few units wide) and sit at
particular frequencies that recur from window to window: 110, 171, 280,
125, 218.

**Those frequencies are the largest gaps between zeta's zeros.**

| window c | W = 2 pi c | state's level | its band | zeros on either side | gap | in mean spacings |
|---|---|---|---|---|---|---|
| 7.2 | 45 | 0.737 | largest mode at 63.7 | 60.83, 65.11 | 4.28 | 1.57 |
| 13.2 | 83 | 0.463 | 107 to 110 | 107.17, 111.03 | 3.86 | 1.75 |
| 16.2 | 102 | 0.398 | 108 to 117 | 107.17, 111.03 | 3.86 | 1.75 |
| 25.2 | 158 | 0.390 | 171 | 169.91, 173.41 | 3.50 | 1.84 |
| 30.2 | 190 | 0.427 | 277 to 284 | 279.23, 282.47 | 3.24 | |
| 20.2 | 127 | 0.758 | 125 to 136 | 124.26, 127.52 | 3.26 | 1.56 |

107.17 to 111.03 is the largest gap between consecutive zeros anywhere
above 66; 169.91 to 173.41 is the largest above 111. So a low ordinary
state is a wave packet parked in an unusually wide gap between zeros:
narrow in frequency, hence long and smooth in the window, hence small at
the window's ends, hence flat. All four properties of steps 2 and 3
follow from that one.

**What this does to the 2/lambda rule.**
- Gaps below W are used up by the tiny states (below W every zero can be
  avoided). So at window c the bottom of the ordinary spectrum is set by
  the best gap above 2 pi c, and it changes hands when W passes a gap:
  the 111 gap serves from about c = 11 to 17.7, the 173 gap from there
  to 27.6, and so on.
- Lowest ordinary level times sqrt(c)/2 at c = 5.2, 7.2, 9.2, 11.2, 13.2,
  16.2, 20.2, 25.2, 30.2: 0.83, 0.86, 0.98, 0.89, 0.84, 0.80, 0.98, 0.98,
  1.07. So 2/lambda holds here to 20%, as Andrews-Connes-La Gatta report
  ("about 2/lambda with fluctuations of about 20%"), and the fluctuations
  are now accounted for: they are which gap is in use.
- But the mechanism does not give 1/lambda. The mean spacing of zeros at
  W is 2 pi / ln c and the window's resolution is 2 pi / L = 2 pi / ln c:
  they shrink together. If the best available gap stays near 1.8 mean
  spacings, (gap) x (window length) stays near 1.8 x 2 pi and the level
  should level off or fall only as slowly as extreme gaps grow. Our
  levels from c = 13 to 30 are 0.46, 0.40, 0.44, 0.39, 0.39: flat, not
  falling like 1/sqrt(c) (which would give 0.55 down to 0.36). Their
  quoted values 0.195 at lambda^2 = 100 and 0.19 to 0.20 at 150 also look
  flat. So 2/lambda may be a fit over a range, and the true law one about
  extreme gaps between zeros. Not settled: we stop at c = 30, N = 300
  limits the frequencies we can see to a few hundred, and we have not
  removed the pole direction as they do.

**Status.** Verified: the identification of the low ordinary states with
the large gaps (six cases above, no counter-example in the nine windows).
Conjecture: the bottom of the ordinary spectrum is governed by the
largest gap between zeros above 2 pi c, measured in units of the window's
resolution. Not found in what we read of Andrews-Connes-La Gatta, who
list "explaining mu about 2/lambda" as open.

**Physics words.** Zeros are absorption lines. A wide gap between two
lines is a transmission window. The lowest excitation above the gas is a
wave packet sitting in the widest transmission window it can find; its
energy is what leaks onto the neighbouring lines.

## Step 5: their range, c = 40 to 150

`flat_track.py` at HP_PREC=128 (the entries are only needed to double
precision; c = 30.2 gives the same eigenvalues as at 420 bits), N = 400,
450, 650, 900, 1150 for c = 40.2, 50.2, 70.2, 100.2, 150.2, THR = 0.03 so
that levels below 0.2 are not mistaken for tiny ones. Zeros up to 2300
from the Riemann-Siegel formula in doubles (`zeta_zeros_to_2300.json`,
1796 found, 1795.6 expected).

| c | W = 2 pi c | lowest ordinary levels (odd half) | where they sit | gap, rank among gaps above W | 2/sqrt(c) |
|---|---|---|---|---|---|
| 40.2 | 253 | 0.329, 0.369, 0.431 | 415.5-418.4, 279.2-282.5, 370.1-373.1 | 2.93 (3rd), 3.24 (1st), 3.01 (2nd) | 0.315 |
| 50.2 | 315 | 0.269, 0.336 | 415.5-418.4, 370.1-373.1 | 2.93 (2nd), 3.01 (1st) | 0.282 |
| 70.2 | 441 | 0.236, 0.248, 0.286 | 415.5-418.4 (just below W), 650.7-653.6, 525.1-527.9 | 2.93, 2.98 (1st), 2.83 (3rd) | 0.239 |
| 100.2 | 630 | **0.1947**, 0.247, 0.266 | 650.7-653.6, 541.8-544.3, 724.6-727.1 | 2.98 (1st) | 0.200 |
| 150.2 | 944 | **0.202**, 0.212, 0.253 | 1068.0-1070.5, 1267.6-1270.1, 1177.1-1179.7 | 2.55 (3rd), 2.55 (2nd), 2.60 (1st) | 0.163 |

- **Agreement with Andrews-Connes-La Gatta's published values.** They
  give mu^odd = 0.194 at lambda^2 = 100 and 0.19 to 0.20 at 150 (for the
  form with the prolate vectors and the pole direction removed). Our raw
  odd levels at c = 100.2 and 150.2: 0.1947 and 0.202.
- **At c = 100 that number is one state in one gap**: a wave packet at
  frequency 652, between the zeros 650.67 and 653.65, the widest gap
  between consecutive zeros above 630 (2.98, about 2.2 mean spacings),
  with 0.3% of its weight in the outer 2% of the window.
- **2/lambda fails at 150** (0.163 predicted, 0.20 found by us and by
  them). The level fell from 0.39 at c = 30 to 0.195 at c = 100 because
  the exceptional gap at 652 came into use (from about c = 66, when W
  passed the gap at 417) and improved as the window lengthened; it is
  used up at c = 104, when W passes it.
- The states are in the top-ranked gaps, but level order is not exactly
  gap order: neighbours matter too. At c = 150.2 our basis only reaches
  frequency 1442, so wider gaps higher up could give a lower level; ours
  is an upper bound there.
- A gap just below W can still host an ordinary state (417 at c = 70.2):
  the tiny states do not use every gap the moment W passes it.
