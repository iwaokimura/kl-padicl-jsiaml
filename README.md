# kl-padicl-jsiaml

SageMath implementation of the Kubota-Leopoldt *p*-adic *L*-function, companion code for:

> Iwao Kimura and Kyota Umeda,
> "Numerical Evaluation of the Kubota-Leopoldt *p*-adic *L*-Function at Positive Integers",
> *JSIAM Letters* (submitted).
> Research activity group: Algorithmic Number Theory and Its Applications.

## Overview

This repository provides three algorithms for computing values of the Kubota-Leopoldt *p*-adic *L*-function L_p(r, χ_F) associated with the Kronecker character χ_F of a quadratic field F = ℚ(√d):

| Algorithm | Method | Best suited for |
|-----------|--------|-----------------|
| **KL** | Kubota-Leopoldt limit formula (partial sum of generalized Bernoulli number definition) | Any positive integer r |
| **p-MT** | p-adic Mellin transform (Riemann-sum approximation of the p-adic integral) | Any integer s (including s ≤ 0) |
| **Diamond** | Diamond's formula using Stirling numbers of the first kind; specialized form for r = 2 | Positive integer r, fastest for r = 2 |

The scripts in `examples/` reproduce the numerical results in the paper.

## Requirements

- [SageMath](https://www.sagemath.org/) (tested with SageMath 10.x)

No additional packages are needed beyond a standard SageMath installation.

## Usage

Run all example scripts from the **repository root**:

```bash
sage examples/table1_lfunction_values.sage
sage examples/table2_execution_times.sage
sage examples/syntomic_regulator.sage
```

### `examples/table1_lfunction_values.sage`

Reproduces **Table 1** of the paper: L_3(2, χ_F) for real quadratic fields
F = ℚ(√d), d ∈ {2, 5, 7, 11}, computed with all three algorithms.
All three algorithms produce the same result.

Expected output (truncated):
```
d=2:  2 + 3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^8 + 3^9 + O(3^10)
d=5:  2 + 2*3^2 + 2*3^3 + 3^6 + 2*3^9 + O(3^10)
d=7:  1 + 2*3 + 3^3 + 2*3^4 + 3^5 + 3^6 + 3^8 + 3^9 + O(3^10)
d=11: 1 + 2*3 + 3^2 + 3^3 + 3^4 + 3^5 + 3^6 + 2*3^7 + 2*3^8 + 3^9 + O(3^10)
```

### `examples/table2_execution_times.sage`

Reproduces **Table 2** of the paper: mean execution time (seconds) of each algorithm
for d ∈ {2, 5, 7, 11}, using 10 repeated runs.
Algorithm Diamond is fastest; Algorithm p-MT is roughly 10²–10³× slower.

**Note on truncation levels**: Both scripts use N = 5 (KL), N = 10 (p-MT), N = 8 (Diamond),
which are the minimum values giving 10-digit accuracy for the implementations here.
The p-MT computations with N = 10 are slow (∼5–30 s per d value, ∼10 min total for
`table2`); set `SKIP_PMT = True` at the top of either script to skip them.

### `examples/syntomic_regulator.sage`

Reproduces the **syntomic regulator table**: R₂ˢʸⁿ(ℚ(√d))/√(−1) for imaginary quadratic fields d ∈ {−1, −2, −5, −10, −11}, computed via

    R₂ˢʸⁿ(F)/√(−1)  =  L_p(2, χ_F·ω⁻¹) · (1 − χ_F(p)·p⁻²)⁻¹ · |d_F|^(3/2)

Expected output (p = 3, 10 digits):
```
d=-1:  3^2 + 2*3^3 + 3^4 + 2*3^5 + 3^7 + 2*3^10 + 2*3^11 + O(3^12)
d=-2:  2*3^3 + 3^4 + 3^5 + 3^7 + 3^8 + 2*3^11 + 3^12 + O(3^13)
d=-5:  2*3^3 + 3^5 + 3^6 + 2*3^7 + 3^8 + 3^9 + O(3^13)
d=-10: 3^2 + 3^3 + 2*3^4 + 3^6 + 3^7 + 2*3^8 + 3^9 + 2*3^10 + 3^11 + O(3^12)
d=-11: 3^3 + 3^5 + 3^6 + 2*3^7 + 3^9 + 3^10 + 2*3^11 + 2*3^12 + O(3^13)
```

## Using the library directly

The core class `KLpAdicL` can be used interactively in SageMath:

```python
sage: load('KLpAdicL.sage')
sage: kl = KLpAdicL(p=3, N=20)          # N = working precision

# Algorithm KL
sage: kl.lfunc_kl_original(2, 2, upto=1)

# Algorithm p-MT
sage: kl.lfunc_padic_integral(2, 2, upto=2)

# Algorithm Diamond (general r)
sage: kl.lfunc_diamond(2, 2, upto=8)

# Algorithm Diamond (specialized r=2, faster)
sage: kl.lfunc_diamond_2_chi(2, upto=8)

# Unified interface (dispatches automatically)
sage: kl.kl_padic_l(2, (3, 2))
```

## License

This project is licensed under the GNU General Public License v3.0 — see [LICENSE](LICENSE) for details.
