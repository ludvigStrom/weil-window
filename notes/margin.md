# Leaving out a prime: the edge response and the margin

Notation: window [1/λ, λ] in the scale variable, c = λ², length
L = ln c in the logarithm of scale, half-length a = L/2. The form is
split into an even and an odd half. ε is an eigenvalue. Everything is
computed in ball arithmetic (python-flint, 420 bits) in a basis of N
sines or cosines; N is given with each result.

Connes and Consani (arXiv:2106.01715, sections 2.2 to 2.4) show in
figures that the form fails without each of the prime powers 2, 3, 4, 5,
7 once the window passes it. These notes ask how soon, and what sets it.

## An exact criterion

For ln n < L < 2 ln n the term of the prime power n only joins two
strips of width w = L - ln n at the two ends of the window. On odd
functions it is minus a reflection inside the strip. With
λ_n = ln p / √n the form without n fails exactly when

    2 λ_n ρ(w) = 1,

where ρ(w) is the top eigenvalue of (strip overlap) relative to the odd
half of the full form: how strongly the two ends of the window respond
to each other. `edge_response.py N c w` computes ρ; the strip matrix is
in closed form and agrees with the withheld term to 1e-18. At the
measured margins of n = 2, 3, 5 the left side is 0.99, 1.02, 1.01.

## Measured margins (odd half, which fails first)

| n | w* | w* n² | basis sizes used | model below |
|---|---|---|---|---|
| 2 | 0.049 | 0.20 | up to 340 | |
| 3 | 0.0163 | 0.147 | up to 480 | 0.0211 |
| 4 | 0.0118 | 0.189 | up to 340 | |
| 5 | 0.0051 | 0.128 | up to 480 | 0.0055 |
| 7 | 0.0025 ± 0.0003 | 0.123 | 240, 480, 720 | 0.0026 |
| 11 | 0.0011 ± 0.0001 | 0.13 | 480, 720, 960 | 0.0010 (predicted first) |

Extrapolated in N; convergence is slow because eigenvectors vanish at
the ends only like 1/√(ln(1/δ)). The n = 11 strip is poorly resolved
even at N = 960. The even half fails later (c = 2.271, 3.140, 4.123,
5.077, 7.058, 8.072).

## What ρ is made of

1. **A floor of about 0.2** from the ordinary states, nearly independent
   of w. Measured up to c = 7, not explained.
2. **Rates.** For each tiny odd eigenvalue, β²/ε, where β is the
   coefficient of 1/√(ln(1/δ)) at the end of the window. By Chekmachev's
   Hadamard formula dε/da = -2β², this is half the rate at which the
   eigenvalue falls as the window grows. From the prolate picture in WKB
   form the rate is 4π²c / K(√q) (K the complete elliptic integral,
   q the prolate eigenvalue over (2πc)²), and the rates sum to π²c².
   Measured sums: 0.85, 0.89, 0.93 of π²c² at c = 7, 11, 13.
3. **Edge shape.** Near the end, at depth t in the logarithm of scale,
   a tiny eigenstate looks like
   β I₀(C√t) / √(ln(1/(κt))), κ = 2π e^γ,
   with C = √(2((2πc)² - χ_n)) from the prolate eigenvalue χ_n. The
   logarithmic factor comes from the symbol ln(|ξ|/2π) on a half-line;
   the Bessel factor is a prolate function at its own endpoint.

Together (`margin_predict.py`):

    2 λ_n [ floor + Σ_j (β_j²/ε_j) ∫₀^w P_j(t) P_j(w - t) dt ] = 1,
    P_j(t) = I₀(C_j √t) / √(ln(1/(κt))).

This gives the last column of the table. For large n the sum closes
(`margin_asymptotic.py`): with μ = w* n²,

    2 λ_n [ floor + (I₀(4π√μ) - 1) / (4 ℓ) ] = 1,   ℓ ≈ ln(1/(κ w*)),

so w* ~ (ln 4√n)² / (16π² n²): the "1/8" in w* n² is a plateau of a
slowly rising function.

## What happens after the margin

Without the term of 7 (odd half, N = 80, which does not resolve the
margin itself): the lowest eigenvalue is -6e-25 at c = 7.03, -2e-18 at
7.04, -7e-10 at 7.06, -8e-5 at 7.10, -0.042 at 7.20, -0.107 at 7.30.
The failure does not stay small.

## Status

- The criterion is exact. The margins are measurements with the stated
  uncertainty. The formula is a model, within 10 to 30% where checked,
  and assumes the form with the earlier primes is positive and has the
  prolate structure. It says nothing about why the margin stays positive.
- Prior work: the prolate structure and the ladder of tiny eigenvalues
  (Connes-Consani; Andrews-Connes-La Gatta, who also give the uniform
  elliptic-integral law and the endpoint value of the prolate
  functions); the Hadamard formula (Chekmachev). Not found elsewhere:
  the criterion, the margins, the model.
- Withdrawn along the way: a "new" identity m′k′ = m² for an associated
  string (it is Chekmachev's Pohozaev identity); a constant correction
  factor of 4 to 5; exponential growth of the deep states from the edge
  (it is the Bessel shape); w* n² = 1/8 exactly.
- A numerical trap: building sin(mθ), cos(mθ) by a rotation recurrence in
  ball arithmetic loses about half a bit per step; above N of about 750
  at 420 bits this silently corrupted the form. `hp_block.py` now
  restarts the recurrence every 32 steps.
