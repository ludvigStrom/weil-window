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
