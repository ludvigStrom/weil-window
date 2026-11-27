# Weil Window

Why this exists

I'm not a mathematician. I like poking at things I don't understand, and like many people I have tried my hand at the Riemann hypothesis. When I heard that Alpöge had announced a proof, I felt some relief, but I couldn't keep my hands off it. I took the argument apart looking for the physics I expected to be hiding in it, and found none. Everything that interested me was missing.

So I turned to Connes's recent papers, where the zeros come from a concrete object: the Weil form on a finite window. I noticed something small, and possibly new, about the odd half of that form, and wanted to push on.

What I wanted was a real physical system whose behaviour is the zeta zeros. I started with cold fermionic gases, because they show the same fine tuning the zeros must have, and soon had a picture of a supercold gas in a box. I began writing a dictionary between Connes's language and the physical language I think in. The background energy in the box read as vacuum energy. The geometry pointed at a horizon, the kind a black hole or an expanding universe has. From there each term in the formula found a place as a physical object: the primes as mirrors, the scaling as a clock, the lowest state as a vacuum.

This repository is what came out: six numerical findings about the Weil form, the scripts and data to check them, and the physical reading that led me to them. None of it is a proof. The picture is a guide, and the numbers are the evidence. This is my attempt to contribute something to the problem.


Numerical experiments on Weil's quadratic form restricted to a finite
window [1/λ, λ] (write c = λ²), the form whose positivity for every
window is equivalent to the Riemann hypothesis.

**What this is.** Observations about a known construction (Connes,
Consani, Moscovici and co-authors), each with its evidence and its
limits. **What it is not.** A proof of anything, or progress towards one
beyond describing more precisely where the difficulty sits. Everything
is computed on finite windows; the hypothesis is the statement about
every window, and nothing here bridges that.

I am not a mathematician. The work was done in conversation with an AI
assistant (Claude), which wrote the code and most of the text; I chose
the questions. Several explanations offered along the way were wrong and
are listed in the notes. The literature was checked lightly, so
"possibly new" means only that I did not find it.

- Plain-words account: [`notes/explain_like_feynman.md`](notes/explain_like_feynman.md)
  (or the simpler [`notes/explain_it_like_im_5.md`](notes/explain_it_like_im_5.md))
- Working log with every number and every dead end: [`notes/questions.md`](notes/questions.md)
- What the picture says the system is: [`notes/clock_and_form.md`](notes/clock_and_form.md)

## Six findings

Write ε for the lowest eigenvalue of a half of the form and call its
eigenvector the vacuum. "Ticks" are the real zeros of the vacuum's
Fourier transform, which approximate the zeros of zeta (Connes, Consani,
Moscovici, arXiv:2511.22755).

**1. The bottom of the ordinary spectrum sits in gaps between zeros.**
The lowest non-minuscule eigenstates are wave packets in wide gaps
between consecutive zeros near 2πc. Since μ^odd of Andrews, Connes and
La Gatta (arXiv:2611.13129) is at most the c-th odd eigenvalue of the raw
form, which needs only a list of zeros (Groskin, arXiv:2607.02828):

| c | c-th odd eigenvalue | their μ^odd | 2/λ | λ·μ^odd at most |
|---|---|---|---|---|
| 50 | 0.269 | 0.270 | 0.283 | |
| 100 | 0.1944 | 0.194 | 0.200 | |
| 250 | 0.097 | | 0.126 | 1.53 |
| 500 | 0.057 | | 0.089 | 1.28 |

Prediction: λ·μ_λ leaves their observed range 1.6 to 2.5 beyond c = 250.
*Limits:* an upper bound beyond c = 150; their values were read through
a summarising tool and should be checked against their Table 7.

**2. The dilation is the only soft way to move the primes, and ε is its
slack.** Move a prime p to p(1+δ): the lowest eigenvalue is −k|δ|, on
both sides. Over all joint moves a blind search finds none that keeps
positivity and exactly one that is ten times softer than the rest:
δ_p ∝ ln p, i.e. p → p^(1+η) (overlap with ln p 0.9989 at c = 14.2,
0.9998 at 28.2). Its cost is archimedean. The zoom-out at which
positivity breaks equals the smallest odd eigenvalue over a rate close
to 1 (7.333e-4 predicted, 7.325e-4 measured at c = 2.6; same to 0.1% at
3.3 and 4.2). *Limits:* first order; primes already in the window.

**3. One offset law and one kernel for the ticks.** t_n = γ_n + K_n·ε,
with K_n tending to a classical number, about 1/(γ_n Ξ'(γ_n)²) (first
zero: measured/limit = 1 − 6.0/c at four windows). A prime power at
lag ln q changes each tick by K_n times the change of ε, which is the
vacuum's autocorrelation at that lag (0.998 of the response at the edge
of the window, over a hundred decades). Near primes follow from
Riemann's Φ alone to 3–8%. *Limits:* first order; a sketch plus tests.

**4. ε against the window: a kick per prime and a note per prime.** Over
322 windows from c = 4 to 24, ε falls steadily; each prime entering
raises the rate of fall by about a decade per unit (5, 7, 11, 13, 17, 19,
23: yes; 6, 10, 12, 14, 15, 18, 20, 21, 22: no); and there are sharp
tones at q − 1/q cycles per unit of c for a mirror at q (2 → 1.50,
3 → 2.65, 4 → 3.74), with pair tones at ratios (3/2, 5/2, 5/3). The
zeros leave no trace in this curve. *Limits:* tone values sit 0.2–0.7%
below the formula; two faint lines unidentified.

**5. The same in two other number systems.** The form for the Dirichlet
characters mod 4 (no prime 2) and mod 3 (no prime 3) is positive, has
the ladder of tiny eigenvalues, and its ticks are the zeros of the
L-function to 27 and 37 digits. Kicks follow the sign of χ(p). The tone
at 1.5 is present without the prime 3 and absent without the prime 2;
the tone at 2.65 the other way round. *Limits:* two characters.

**6. A window-only estimate of a tick's distance from its limit.** With
B_n = t_n P_n ε / v̂'(t_n)², P_n = Σ v_m²/(t_n² − ω_m²), everything
computed on the window: (t_n − γ_n)/B_n lies between 1 and
1 + 0.0035·(γ_n N/c)² in all 64 cases tried (three systems, four zeros,
sixteen windows). A conjecture to aim for. *Limits:* the identity behind
it assumes the zeros are real; it describes the convergence if the
hypothesis holds, it does not prove it.

## Reproducing it

Python 3.10 or later with `python-flint`, `mpmath`, `numpy`, `scipy`.
Scripts read and write in the current folder, so run them from `data/`.

```
pip install python-flint mpmath numpy scipy
cd data

# 1  levels from the zeros alone
python3 ../scripts/mu_from_zeros.py 2.0 100.2 148.2
python3 ../scripts/mu_big.py 0 2 3.2 500.2

# 2  moving the primes
HP_PREC=128 python3 ../scripts/prime_move.py 100 14.2 7 1e-3 -1e-3 1e-2
HP_PREC=128 python3 ../scripts/rigid.py 100 14.2
HP_PREC=128 python3 ../scripts/zoom_slack.py 60 3.3

# 3  ticks, offset law, kernel
HP_PREC=300 python3 ../scripts/clock_ticks.py 50 8.2 60
HP_PREC=600 HP_DIGITS=80 python3 ../scripts/sens_fast.py 70 14.2
python3 ../scripts/near_primes.py
python3 ../scripts/k_multiplier.py

# 4  the level curve (stored in level_curve.json) and its tones
python3 ../scripts/bumps.py

# 5  another number system
HP_PREC=420 HP_DIGITS=60 DPS=70 python3 ../scripts/lticks.py chi4 90 24.2 22

# 6  the window-only estimate
HP_PREC=560 HP_DIGITS=100 python3 ../scripts/tick_bound.py 90 8.2
```

`hp_block.py` builds the form in ball arithmetic; `HP_PREC` is the
working precision in bits and `HP_DIGITS` the accuracy asked of the
quadrature. The smallest eigenvalue at window c is about 10^(−5c), so the
precision has to grow with the window.

## Prior work this rests on

- A. Connes, C. Consani, *Spectral triples and ζ-cycles*, arXiv:2106.01715
- A. Connes, C. Consani, H. Moscovici, *Zeta spectral triples*, arXiv:2511.22755
- A. Andrews, A. Connes, G. La Gatta, arXiv:2611.13129
- A. Connes, H. Moscovici, arXiv:2112.05500
- A. Groskin, arXiv:2607.02828, arXiv:2605.20224
- E. Bombieri, *Remarks on Weil's quadratic functional in the theory of prime numbers, I* (2000)
- M. Suzuki, arXiv:2610.41115, arXiv:2606.09096; J.-F. Burnol, arXiv:math/9809119
- D. Slepian, H. Pollak; W. Fuchs: prolate functions and their asymptotics

## Licence and contact

MIT (see `LICENSE`). Ludvig Ström. Corrections are very welcome.
