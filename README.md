# weil-window

Numerical experiments on Weil's quadratic form restricted to a finite
window [1/λ, λ] (write c = λ²), the form whose positivity for all windows
is equivalent to the Riemann hypothesis.

This is a hobby project. I am not a mathematician. The work was done in
conversation with an AI assistant (Claude), which wrote the code and most
of the text; I chose the questions and checked what I could. Nothing here
is a proof of anything, and mistakes are possible.

Start with [`notes/explain_it_like_im_5.md`](notes/explain_it_like_im_5.md)
for the whole story in plain words. The working log, with every number
and every dead end, is [`notes/questions.md`](notes/questions.md).

There are two observations.

## 1. The bottom of the ordinary spectrum sits in gaps between zeros

The form has about 2c extremely small eigenvalues (Connes and Consani)
and then an ordinary part. Andrews, Connes and La Gatta
([arXiv:2611.13129](https://arxiv.org/abs/2611.13129)) find that the
bottom of the ordinary part, μ_λ, is close to 2/λ, with fluctuations of
about 20%, for λ² from 5 to 150, and list explaining this as open.

What the computations here show:

- **The lowest ordinary eigenstates are wave packets sitting in wide gaps
  between consecutive zeta zeros near 2πc.** At c = 100 it is one state at
  frequency 652, between the zeros 650.67 and 653.65. As the window grows
  the level follows one gap for a while, then hands over to another.
- **These levels need only a list of zeros.** For an odd test function
  the form is twice the sum over zeros of the squared sine transform at
  the zero, so in a sine basis it is 2AAᵀ with A in closed form. (That the
  truncated form is an exact sum over zeros is in A. Groskin,
  [arXiv:2607.02828](https://arxiv.org/abs/2607.02828); here it is only
  used.) It reproduces a ball-arithmetic computation of the form to three
  or four digits.
- **A prediction.** Their μ^odd is defined by removing λ² − 1 directions,
  so it is at most the λ²-th odd eigenvalue of the raw form:

| λ² | λ²-th odd eigenvalue, from the zeros | their μ^odd | 2/λ | λ·μ^odd at most |
|---|---|---|---|---|
| 28 | 0.374 | 0.361 | 0.378 | |
| 50 | 0.269 | 0.270 | 0.283 | |
| 100 | 0.1944 | 0.194 | 0.200 | |
| 150 | about 0.20 | 0.194 | 0.163 | |
| 250 | 0.097 | | 0.126 | 1.53 |
| 300 | 0.079 | | 0.115 | 1.37 |
| 500 | 0.057 | | 0.089 | 1.28 |

  Their reported range is 1.6 ≤ λμ ≤ 2.5. If these numbers are right,
  λμ leaves that range beyond λ² = 250, and 2/λ is a description of the
  range computed so far, not the long-run rule.

```
cd data
python3 ../scripts/level_from_zeros.py 100.2 900        # from zeros vs the stored ball-arithmetic run
python3 ../scripts/mu_from_zeros.py 2.0 100.2 148.2     # lowest ordinary level and its gap (zeros to 6000)
python3 ../scripts/mu_big.py 0 2 3.2 500.2              # large windows (zeros to 21000)
```

Caveats: beyond 150 this is an upper bound and is not shown to stay
tight; odd half only; windows at c + 0.2; zeros from the Riemann–Siegel
formula in double precision; their numbers were read through a
summarising tool and should be checked against the paper.

## 2. Moving the primes

What happens to the form when the primes are not where they are?

- **A tent.** Move one prime p to p(1+δ): the lowest eigenvalue is
  −k·|δ|, negative on both sides, linear from the first step. Positivity
  holds at the true position only. Exactly, k = (ln p/√p) × the extreme
  eigenvalue, on the near-null states, of the derivative of the shift
  operator at distance ln p.
- **One soft direction.** A blind search over all joint moves of the
  primes finds no direction that keeps positivity, and exactly one that is
  ten times softer than the rest: every prime moved in proportion to
  ln p, that is p → p^(1+ε), the dilation of the log axis.
- **The background sets the scale.** The cost of that zoom comes almost
  entirely from the archimedean and pole terms, which do not stretch.
  Zooming out sinks all near-null levels together; zooming in lifts them
  all except one state, which is the pole vector (97%).
- **The tiny eigenvalues are zoom slack.** The smallest odd eigenvalue,
  divided by a rate close to 1, is the zoom-out at which positivity
  breaks, checked to 0.1% at c = 2.6, 3.3, 4.2. So a window up to c fixes
  the overall scale of the primes to about e^(−4πc).
- **Curves are different.** The same test on the window form of a curve
  over a finite field (where positivity is Weil's theorem) finds 2g − 1
  exactly free directions and an open set of allowed ones. The primes,
  within the moves tried, have none.

```
cd data
HP_PREC=128 python3 ../scripts/prime_move.py 100 14.2 7 1e-3 -1e-3 1e-2   # the tent
HP_PREC=128 python3 ../scripts/rigid.py 100 14.2                          # softest direction, softness spectrum
HP_PREC=128 python3 ../scripts/zoom_pole.py 100 14.2 0 -1e-2 1e-2         # zoom in / out, with and without the pole
HP_PREC=128 python3 ../scripts/zoom_slack.py 60 3.3                       # smallest eigenvalue against critical zoom
python3 ../scripts/ff_rigid.py 3 7 0.3041 0.6959                          # the curve toy
```

Caveats: windows 14.2 and 28.2 (2.6 to 4.2 for the slack test), one
curve; first-order tests, so small moves only; only primes already in the
window are moved. That wrong primes break positivity is expected by
everyone; a light literature search found nobody measuring it this way,
which is not the same as nobody having done it.

Two attempts to reduce the near-null states to a small hopping model (by
place, then by frequency) failed by tests fixed beforehand; see the log.

## Requirements

Python 3.10 or later with `python-flint`, `mpmath`, `numpy`, `scipy`.
Scripts read and write in the current folder, so run them from `data/`.
`rigid.py` writes a cache (`rigid_P_*.npz`, not in the repository) that
`rigid_weights.py`, `stiff_profile.py`, `chain_model.py` and
`pitch_model.py` reuse.

## What is in the repository

| Path | What |
|---|---|
| `scripts/hp_block.py` | Weil's form on a window, even and odd halves, in ball arithmetic |
| `scripts/level_from_zeros.py`, `mu_from_zeros.py`, `mu_big.py`, `zeros_rs.py` | Levels from the zeros alone; zeta zeros by Riemann–Siegel |
| `scripts/flat_track.py`, `gap_table.py` | The gap observation from the form itself |
| `scripts/prime_move.py`, `rigid.py`, `rigid_weights.py`, `soft.py`, `neg_count.py` | Moving primes: tent, rigidity, counts of negative levels |
| `scripts/zoom_decomp.py`, `zoom_pole.py`, `zoom_slack.py` | The zoom: where its cost comes from, the pole state, slack |
| `scripts/stiff_profile.py` | The exact stiffness rule and its profile across the window |
| `scripts/ff_toy.py`, `ff_rigid.py`, `diver_floor.py` | Curves over finite fields; leaving primes out |
| `scripts/chain_model.py`, `pitch_model.py` | The two failed small models |
| `scripts/ladder_track.py`, `edge_*.py`, `deep_profile.py`, `margin_*.py` | Earlier work: births of tiny levels, leaving out one prime, the margin |
| `scripts/horizon_*.py`, `horn_*.py` | Earlier work: checks of a physical reading |
| `data/` | Outputs quoted in the notes; zeta zeros to 2300, 6000, 21000 |
| `notes/explain_it_like_im_5.md` | The whole project in plain words |
| `notes/questions.md` | Working log for both observations |
| `notes/ladder.md`, `margin.md`, `picture.md`, `function_field.md` | Earlier notes |

## Prior work this rests on

- A. Connes, C. Consani, *Spectral triples and ζ-cycles*, arXiv:2106.01715
- A. Andrews, A. Connes, G. La Gatta, arXiv:2611.13129
- A. Connes, H. Moscovici, arXiv:2112.05500; A. Connes, C. Consani, H. Moscovici, arXiv:2511.22755
- A. Groskin, arXiv:2607.02828 and arXiv:2605.20224
- E. Bombieri, *Remarks on Weil's quadratic functional in the theory of prime numbers, I* (2000)
- M. Suzuki, arXiv:2610.41115 and arXiv:2606.09096; J.-F. Burnol, arXiv:math/9809119
- D. Slepian, H. Pollak; W. Fuchs: prolate functions and their asymptotics

## Contact

Ludvig Ström. Corrections are very welcome.
