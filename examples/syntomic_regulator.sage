"""
Reproduce the syntomic regulator table of Kimura & Umeda, JSIAM Letters (submitted).

Computes R_2^syn(F) / sqrt(-1) for imaginary quadratic fields F = Q(sqrt(d)),
d in {-1, -2, -5, -10, -11}, using the formula

    R_2^syn(F) / sqrt(-1)  =  L_p(2, chi_F * omega^{-1})
                               * (1 - chi_F(p) * p^{-2})^{-1}
                               * |d_F|^{3/2}

where d_F is the discriminant of F, chi_F is the Kronecker character (d_F / ·),
omega is the Teichmüller character mod p, and L_p denotes the Kubota-Leopoldt
p-adic L-function.

Reference:
  Kolster, M., Nguyen Quang Do, T. and Fleckinger, V., "Twisted S-units,
  p-adic class number formulas, and the Lichtenbaum conjectures",
  Duke Math. J. 84 (1996).

Usage (run from the repository root):
    sage examples/syntomic_regulator.sage
"""

# ── Helper functions ─────────────────────────────────────────────────────────

def _chid_omegap(d, p, r, pre):
    """Return chi_d * omega_p^r as a Dirichlet character over pAdicField(p, pre)."""
    Kp = pAdicField(p, pre)
    DGp = DirichletGroup(p, Kp)
    e = abs(4*d) if d % 4 in [2, 3] else abs(d)
    DGep = DirichletGroup(e*p, Kp)
    chid = kronecker_character(d)
    omegap = DGp.gen()
    return DGep(chid) * DGep(omegap)^r


def lfunc_diamond_2_chiomega(p, d, N, pre):
    """Compute L_p(2, chi_d * omega^{-1}) via Diamond's formula specialized to r=2.

    Uses the series

        L_p(2, chi * omega^{-1}) = sum_{a, gcd(a,p)=1}^{p*f}
            chi(a) * sum_{n>=0} n! / (n+1) * (p*f)^{n-1}
                     / (a * (a + p*f) * ... * (a + n*p*f))

    where f is the conductor of chi = kronecker_character(d) and the omega^{-1}
    twist is absorbed (chi itself, without the omega factor, appears in the sum).

    Parameters
    ----------
    p : int
        An odd prime.
    d : int
        A squarefree negative integer with gcd(d, p) = 1.
    N : int
        Truncation level of the inner series.
    pre : int
        Working p-adic precision.

    Returns
    -------
    p-adic number
        Approximation of L_p(2, chi_d * omega^{-1}).
    """
    if not is_squarefree(d):
        raise ValueError("d must be squarefree.")
    if d % 4 == 2 or d % 4 == 3:
        d = 4 * d
    K = pAdicField(p, pre)
    chi = kronecker_character(d)
    f = chi.conductor()
    return K(sum(
        chi(a) * sum(
            factorial(n) * (p*f)^(n-1) / ((n+1) * product(a + i*p*f for i in range(n+1)))
            for n in range(N+1)
        )
        for a in range(1, p*f + 1) if gcd(a, p) == 1
    ))


def R_syntomic(p, d, N, pre, dig):
    """Compute R_2^syn(Q(sqrt(d))) / sqrt(-1) for an imaginary quadratic field.

    Parameters
    ----------
    p : int
        An odd prime.
    d : int
        A squarefree negative integer with gcd(d, p) = 1.
    N : int
        Truncation level for the Diamond series (see lfunc_diamond_2_chiomega).
    pre : int
        Working p-adic precision.
    dig : int
        Number of p-adic digits shown in the output.

    Returns
    -------
    p-adic number
        Approximation of R_2^syn(Q(sqrt(d))) / sqrt(-1) at dig-digit precision.
    """
    chid = kronecker_character(d)
    e = abs(4*d) if d % 4 in [2, 3] else abs(d)   # |d_F|, the absolute discriminant
    K = pAdicField(p, pre)
    K_dig = pAdicField(p, dig)
    lfunc = lfunc_diamond_2_chiomega(p, d, N, pre)
    rv = lfunc / (1 - chid(p) * p^(-2)) * K(e)^K(3/2)
    return K_dig(rv)


# ── Parameters ───────────────────────────────────────────────────────────────
p   = 3
N   = 10   # truncation level (sufficient for dig-digit output)
pre = 20   # working p-adic precision
dig = 10   # output digits

d_list = [-1, -2, -5, -10, -11]

# ── Computation ──────────────────────────────────────────────────────────────
print(f"R_2^syn(Q(sqrt(d))) / sqrt(-1)  for imaginary quadratic fields")
print(f"p = {p},  N = {N},  computation precision = {pre},  output digits = {dig}")
print()
print(f"{'d':>5}  {'Result'}")
print("-" * 70)

for d in d_list:
    val = R_syntomic(p, d, N, pre, dig)
    print(f"{d:>5}  {val}")
