# A physical reading of the window form

A dictionary between Weil's form on a finite window and a physical
system, with the numerical checks made of it. It is a reading of known
mathematics, not a result, and it does not explain why the form is
positive.

## The dictionary

| In the Weil form | In the picture | Source of the mathematics |
|---|---|---|
| Window of scales 1/λ to λ | A box; or the visible patch of a universe with a horizon | Connes-Consani |
| Cut-off in scale and in frequency | Positions inside the box; momenta inside a Fermi sea | Connes-Consani; Slepian-Pollak |
| Scaling | The clock | Berry-Keating; Connes |
| Prolate functions | Modes of the gas lying almost entirely inside the box | Slepian-Pollak; Eisler-Peschel for the gas |
| Tiny eigenvalues, about 2c of them | How much each such mode leaks out | Connes-Consani; Andrews-Connes-La Gatta |
| Ends of the window | Horizons, for modes near the surface of the gas | Connes-Moscovici (a remark) |
| Term of the prime power n, weight ln p / √n | A partly reflecting mirror | Weil's explicit formula; cf. Sierra's mirrors |
| ln p and 1/√n | Energy of a "primon"; its thermal factor at inverse temperature ½ | Julia; Spector; Bost-Connes |
| Zeros of zeta | Frequencies the dressed wave cannot carry | Connes |
| Low ordinary eigenstates | Wave packets in the widest gaps between those frequencies | notes/ladder.md |
| Positivity of the form | Stability of the vacuum | not explained |

## Checks

**The count.** A window of length ln c has 2c more modes below frequency
2πc than zeta has zeros of size below 2πc. Those are the tiny states. In
each half one appears per unit of c: at c = 1.96, 2.97, 3.98, 4.98, 6.00,
6.99, 8.01 in the odd half and 2.45, 3.48, 4.46, 5.49, 6.45, 7.50, 8.48
in the even half, whether or not the integer is prime
(`ladder_track.py`).

**The horizon's temperature** (`horizon_temperature.py`). If the ends are
horizons for the surface of the gas, the largest tiny eigenvalues should
be Boltzmann factors in the prolate eigenvalue χ_n:
ln(1/ε) = (π/2W)(W² - χ_n), W = 2πc, n = 6, 10, 14, ... for the odd
half, with nothing fitted. At c = 13.01: measured 24.6, 16.4, 8.3, 2.4
for the four shallowest, law 24.3, 16.0, 8.4, 2.0. The deepest is 126.1
against 108.4: deep states are colder than thermal (e^(-4πc), not
e^(-π²c)). The uniform elliptic-integral law of Andrews-Connes-La Gatta
gives every level to about one unit in the logarithm (127.3 for the
deepest), so this is also a numerical check of the odd half of their
Hypothesis 3.4.

**Echoes** (`horizon_image.py`, `horizon_scan.py`). In the eigen-equation
the term of n places, at depth L - ln n below one end, a copy of the
state's profile at the other end. Measured at c = 6.5 in units of each
state's edge amplitude: 0.08, 0.047, 0.07 at the depths of 5, 4, 3
(weights 0.72, 0.35, 0.63); 0.06 at depth ln 2; 0.023 at the depth of 6,
which has no term of its own (2 then 3); about 0.01 at ln(5/4),
ln(4/3), ln(3/2), ln(5/3), ln(5/2); below 0.003 at five control depths.
Single echoes are 0.11 to 0.14 times the weight, double echoes 0.03 to
0.04 times the product. The echo points are the ends multiplied by
rationals: the identification made in Connes's adele class space.

**The bath.** The mean energy of a gas with one kind of particle per
prime, energy ln p, is Σ Λ(n) n^(-β). At β = ½ this is the weight in the
form. β = ½ is forced: rescaling by n with the norm kept changes the
amplitude by √n. Composites have no term because the form sees the mean
energy; the double echoes multiply as they must for independent
particles. That this gas and the space carrying the zeros are dual is
Bost-Connes and Connes-Consani-Marcolli.

**A gravitational system with the same parts** (`horn_wave.py`).
Hartnoll and Yang (arXiv:2502.02661): near a spacelike singularity the
metric at each point is a particle in half the modular domain; with
Dirichlet walls the states are odd cusp forms, whose L-functions appear
on the critical line in the dilatation basis and are partition functions
of primon-type gases. With Neumann walls the Eisenstein series are
allowed, and standard formulas give for the scattering state of energy ε
the transform ξ(½ + i(t + ε)) ξ(½ + i(t - ε)): zeta. Computed for ε = 2.5
from divisor sums and Bessel functions alone: reflection coefficient
0.804465436047 - 0.593999463136 i, of modulus 1 to 12 digits and equal to
ξ(2iε)/ξ(1 + 2iε); fifteen dark lines at the zeros ± 2.5. This checks a
classical identity (Hecke; Zagier 1981).

**A cut horn does not help** (`horn_window.py`). Lagarias and Suzuki
prove that the horn cut at height T ≥ 1 has all its zeros on the line.
Our window is the same cut with T = √c, but the transform of the wave
over the window is exactly ξξ + (wall term) - (tail), and the wall term
is 0.0077 at the first zero where ξξ is 8e-10. The theorem is about the
wall.

## Horizon, vacuum, dark energy

The Connes-Moscovici metric is ds² = -α dt² + dx²/α with α = 4(λ² - x²):
positive inside the window, horizons at ±λ, scalar curvature +8. That is
the static patch of two-dimensional de Sitter space, the geometry of an
observer surrounded by a cosmological horizon. They call it a black hole
and do not use the words de Sitter; in two dimensions the descriptions
overlap. In that language: the form is positive = the vacuum is stable;
leaving out a prime produces an instability (notes/margin.md); the
bottom of the ordinary spectrum is a zero-point cost, set by the widest
gap between absorption lines (notes/ladder.md).

## Ledger

**Known:** everything in the right-hand column of the dictionary; the
ladder of tiny eigenvalues and its uniform law; the failure of the form
without a small prime (figures in Connes-Consani); zeta from the
Eisenstein channel; causality of an adelic scattering equivalent to RH
(Burnol).

**Not found elsewhere, after a light search:** the two cut-offs read as
the entanglement problem of a Fermi gas; the thermal check; the echoes;
the edge response and the margin; the low ordinary states in the gaps
between zeros; the singularity billiard with Neumann walls read as
zeta's setting.

**Not explained by any of it:** positivity.
