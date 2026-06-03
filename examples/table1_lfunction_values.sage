"""
Reproduce Table 1 of Kimura & Umeda, JSIAM Letters (submitted).

Computes L_3(2, chi_F) for real quadratic fields F = Q(sqrt(d)),
d in {2, 5, 7, 11}, using all three algorithms (KL, p-MT, Diamond).
All three algorithms produce identical results.

Truncation levels used (minimum N for 10-digit accuracy):
  Algorithm KL      : N = 5  (sums f * 3^5 terms, where f = conductor of chi_F)
  Algorithm p-MT    : N = 10 (sums d' * 3^10 terms; takes ~5–30 s per d value)
  Algorithm Diamond : N = 8  (inner series truncated at n = 8)

Note: Algorithm p-MT with N = 10 requires several minutes to finish for
      all four values of d.  To skip it, set SKIP_PMT = True below.

Usage (run from the repository root):
    sage examples/table1_lfunction_values.sage
"""

SKIP_PMT = False   # set True to skip Algorithm p-MT (saves time)

# ── Algorithm implementations ─────────────────────────────────────────────────

def _chid_omegap(d, p, r, pre):
    """Return chi_d * omega_p^r as a Dirichlet character over pAdicField(p, pre)."""
    Kp = pAdicField(p, pre)
    DGp = DirichletGroup(p, Kp)
    e = abs(4*d) if d % 4 in [2, 3] else abs(d)
    DGep = DirichletGroup(e*p, Kp)
    chid = kronecker_character(d)
    omegap = DGp.gen()
    return DGep(chid) * DGep(omegap)^r


def _coprime_to(n):
    """Return the smallest integer >= 2 coprime to n."""
    a = 2
    while gcd(a, n) != 1:
        a += 1
    return a


def lfunc_KL(p, r, d, N, pre, dig):
    """Algorithm KL: N-th partial sum of the Kubota-Leopoldt limit formula."""
    f = kronecker_character(d).conductor()
    chiomega = _chid_omegap(d, p, r-1, pre)
    K = pAdicField(p, pre)
    rv = (1/K((r-1)*f*p^N)) * sum(
        chiomega(a) * K(a)^(1-r)
        for a in range(1, f*p^N + 1) if gcd(a, p) == 1
    )
    return pAdicField(p, dig)(rv)


def lfunc_pMT(p, s, d, N, pre, dig):
    """Algorithm p-MT: N-th Riemann sum of the p-adic Mellin transform."""
    chiomega = _chid_omegap(d, p, -1, pre)
    chi = kronecker_character(d)
    f = chi.conductor()
    d_prime = f // p if f % p == 0 else f
    K = pAdicField(p, pre)
    alpha = _coprime_to(d_prime * p^N)
    S = IntegerModRing(d_prime * p^N)
    ta = K.teichmuller(alpha)
    C = -(ta^(1-s)) / (ta^(1-s) - chi(alpha)*K(alpha)^(1-s))
    inv_alpha = ZZ(S(alpha)^(-1))
    rv = C * sum(
        chiomega(a)
        * (K.teichmuller(a + p^N) / K(a + p^N))^s
        * (K((alpha-1)/2) + alpha/K(d_prime*p^N)
           * (a*inv_alpha - S(a*inv_alpha).lift()))
        for a in range(d_prime * p^N) if gcd(a, p) == 1
    )
    return pAdicField(p, dig)(rv)


def lfunc_Diamond(p, d, N, pre, dig):
    """Algorithm Diamond: N-th partial sum of Diamond's formula for r=2."""
    chiomega = _chid_omegap(d, p, 1, pre)
    f_tilde = chiomega.conductor()
    K = pAdicField(p, pre)
    pf = p * f_tilde
    rv = sum(
        chiomega(a) * sum(
            K(factorial(n)) * K(pf)^(n-1)
            / K((n+1) * product(a + i*pf for i in range(n+1)))
            for n in range(N+1)
        )
        for a in range(1, pf+1) if gcd(a, p) == 1
    )
    return pAdicField(p, dig)(rv)


# ── Parameters ────────────────────────────────────────────────────────────────
p   = 3
r   = 2
pre = 20   # working p-adic precision
dig = 10   # number of 3-adic digits shown

# Minimum truncation levels N for dig-digit accuracy
N_KL  = 5    # valid for all d in {2, 5, 7, 11}
N_pMT = 10   # valid for all d in {2, 5, 7, 11}
N_Dia = 8    # valid for all d in {2, 5, 7, 11}

d_list = [2, 5, 7, 11]

# ── Computation ───────────────────────────────────────────────────────────────
print(f"L_{p}({r}, chi_F)  for real quadratic fields F = Q(sqrt(d))")
print(f"p = {p},  computation precision = {pre},  output digits = {dig}")
if SKIP_PMT:
    print("(Algorithm p-MT skipped; set SKIP_PMT = False to include it)")
print()

col = 54
header = f"{'d':>4}  {'Alg. KL':<{col}}  {'Alg. Diamond':<{col}}"
if not SKIP_PMT:
    header += f"  {'Alg. p-MT':<{col}}"
print(header)
print("-" * len(header))

for d in d_list:
    val_KL  = lfunc_KL(p, r, d, N_KL, pre, dig)
    val_Dia = lfunc_Diamond(p, d, N_Dia, pre, dig)
    if not SKIP_PMT:
        val_pMT = lfunc_pMT(p, r, d, N_pMT, pre, dig)
        print(f"{d:>4}  {str(val_KL):<{col}}  {str(val_Dia):<{col}}  {str(val_pMT):<{col}}")
    else:
        print(f"{d:>4}  {str(val_KL):<{col}}  {str(val_Dia):<{col}}")
