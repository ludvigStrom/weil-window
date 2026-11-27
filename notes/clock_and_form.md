# The clock and the form: what the construct says the system is

26 November 2026. Why we care about Weil positivity at all: we want the
operator behind it. The form is the way in. This note states what our
physical reading says that operator and its setting are, layer by layer,
with the confidence of each layer. Plain words; numbers are in
`questions.md` and `shape_brainstorm.md`.

## In one paragraph

A small expanding universe with one space and one time dimension, seen
by a single observer with a horizon on each side. The operator we are
after is that observer's clock: the generator of the zoom. The form is
the energy of the empty space he sees, and its lowest state fixes how
waves reflect at his horizons. With that reflection rule, the clock's
allowed frequencies are the zeros.

## The layers

| Layer | What it says | Status |
|---|---|---|
| The clock is the zoom | The one soft motion of the primes is p -> p^(1+eps), the stretch of the log ruler. In an expanding universe that stretch is the observer's time. So the operator is the generator of dilations. | Ours as a measurement (blind search, `rigid.py`); Connes's starting point as theory. |
| It expands and has horizons | The window's edges behave like the horizons of two-dimensional de Sitter space (metric -a dt^2 + dx^2/a, a = 4(lambda^2 - x^2), curvature +8). The top of the ladder of tiny levels follows a heat law with no free parameter. | Metric: Connes-Moscovici. Reading as a static patch and the heat-law check: ours. |
| The near-null states are the vacuum | Same structure as empty space seen through a finite patch (prolate functions as entanglement modes). The tiny levels form a geometric ladder, as Boltzmann factors do. | Eisler-Peschel for the free gas; applying it here: ours. |
| The horn is the same object without time | Hyperbolic plane and de Sitter space share the symmetry SL(2,R). On the modular surface a wave sent into the cusp returns with a phase made of zeta. | Classical (Hecke, Zagier, Lax-Phillips); reflection checked numerically (`horn_wave.py`). |
| The primes are images | The horn is a folded space; the primes arise from the folding, as images in mirrors. Each prime touches every near-null state at once (`jam_test.py`), as an image would. | Speculation; the collective action is measured. |
| Projected, not real | Near a cosmological singularity a 3+1 dimensional universe's shape moves as a billiard on this horn. Our system would be the shadow of that. | Hartnoll-Yang; borrowed. |
| The background calibrates the clock | The zoom costs only through the archimedean and pole terms. The tiny levels are the tolerance on the clock rate, closing like e^(-4 pi c). | Ours (`zoom_decomp.py`, `zoom_slack.py`). |

## How the form leads to the operator

Connes, Consani and Moscovici (arXiv:2511.22755) and Connes-van Suijlekom:
take the lowest state of the form on the window; the zeros of its
Fourier transform are real and approximate the zeros of zeta; they are
the spectrum of the scaling operator on the window with a boundary rule
(a rank-one change) built from that state. In our clothes: the vacuum
decides how the horizon reflects, and the clock then ticks at the zeros.

So the form is not the goal. It fixes the boundary.

## The known wall

On the horn the zeros appear as fading echoes (resonances at half the
zeros), not steady tones. Making them steady is what the right boundary
rule has to do, and it is only known window by window.

## The experiment this suggests

Run the recipe ourselves: lowest state, its transform, its zeros. Then
change the tape speed (zoom the primes) and move a single prime, and
watch the ticks. Results in `questions.md`, section "The clock's ticks".

## Added 27 November: what the experiments said

The experiment proposed above was run, and four more after it. Numbers
are in `questions.md` from "The clock's ticks" onwards.

| Finding | In plain words | Status |
|---|---|---|
| The recipe works | With the primes 2, 3, 5, 7 the clock ticks at 13 zeros, the first to 30 digits; with primes up to 23, to 126 digits. | Reproduces Connes-Consani-Moscovici. |
| Tick error is 10^4 times the smallest level | The tiny level is the accuracy of the zeros, as well as the play in the clock rate and the closeness to failure. | Ours, four windows. |
| The clock is the zoom | Zoom the primes by eps and the ticks sit at zero / (1 + eps). | Ours as a measurement. |
| Wrong primes bring back the box | Any move puts extra ticks below the first zero, spaced like the empty window's own modes. Right primes cancel them. | Ours. |
| The vacuum outlives positivity | Moving an outer prime makes the form negative but leaves the vacuum 99.99% intact and the first zero right to 19 digits. Inner primes or a zoom destroy it. | Ours. |
| One kernel | A mirror at any distance q, prime or not, moves a zero in proportion to the intensity of one wave (the sixth of the prolate family) at sqrt(q). Only sqrt(2) is inside its bright part. | Ours; fitted to 0.5 decade over 65. |
| The vacuum is the template | Each zero's response is the vacuum with that zero divided out; for far mirrors the zero's own pitch drops out. | A sketch with a supporting test; a slowly varying factor is unexplained. |

**What this adds to the picture of the operator.**

- The construction is front-loaded: the background and the primes 2 and 3
  make the low zeros; each further prime adds an exponentially thinner
  correction, read off a fixed wave. Good news for a limit existing; not
  a proof that it does.
- The zeros are not added to an empty box by the primes. They are what
  is left when the primes silence the box.
- One number per window, the depth of the vacuum, governs the clock's
  tolerance, the zeros' accuracy, and the scale on which primes fade in.
- The objects to think with are not the individual near-null states
  (their contributions cancel to tens of digits) but one wave, the
  vacuum, and its versions with a zero removed.

**Still the wall.** All of this is on finite windows. Nothing here says
why the vacuum exists at every window.
