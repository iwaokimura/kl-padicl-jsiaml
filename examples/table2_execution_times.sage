"""
Reproduce Table 2 of Kimura & Umeda, JSIAM Letters (submitted).

Measures the execution time of the three algorithms (KL, p-MT, Diamond)
for computing L_3(2, chi_F), F = Q(sqrt(d)), d in {2, 5, 7, 11}.

Truncation levels used (minimum N for 10-digit accuracy):
  Algorithm KL      : N = 5  (f * 3^5 terms; fast, < 0.1 s per d)
  Algorithm p-MT    : N = 10 (d' * 3^10 terms; slow, 5–30 s per d)
  Algorithm Diamond : N = 8  (fast, < 0.1 s per d)

WARNING: This script takes approximately 10–20 minutes to complete,
         mostly due to the p-MT measurements.  To skip p-MT timing,
         set SKIP_PMT = True below.

Usage (run from the repository root):
    sage examples/table2_execution_times.sage
"""

import time as _time

SKIP_PMT = False   # set True to skip Algorithm p-MT timing

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


# ── Timing utility ────────────────────────────────────────────────────────────

def measure(func, args, n_repeat=10):
    """Return (mean_sec, std_sec) over n_repeat runs of func(*args)."""
    times = []
    for _ in range(n_repeat):
        t0 = _time.perf_counter()
        func(*args)
        times.append(_time.perf_counter() - t0)
    mean = sum(times) / n_repeat
    std = (sum((t - mean)^2 for t in times) / n_repeat)^(1/2)
    return mean, std


# ── Parameters ────────────────────────────────────────────────────────────────
p     = 3
r     = 2
pre   = 20   # working p-adic precision
dig   = 10   # output digits
n_rep = 10   # timing repetitions (reduce for p-MT to save time)

# Minimum N for dig-digit accuracy (same values as table1_lfunction_values.sage)
N_KL  = 5
N_pMT = 10
N_Dia = 8

d_list = [2, 5, 7, 11]

# ── Measurement ───────────────────────────────────────────────────────────────
print(f"Execution time comparison: L_{p}({r}, chi_F), F = Q(sqrt(d))")
print(f"p={p}, pre={pre}, dig={dig}")
print(f"N_KL={N_KL}, N_pMT={N_pMT}, N_Diamond={N_Dia}, repetitions={n_rep}")
if SKIP_PMT:
    print("(Algorithm p-MT skipped; set SKIP_PMT = False to include it)")
print()

rows = {}
for d in d_list:
    print(f"d={d}:")
    print(f"  KL (N={N_KL})...", end=" ", flush=True)
    m_KL, s_KL = measure(lfunc_KL, (p, r, d, N_KL, pre, dig), n_rep)
    print(f"{m_KL:.5f} s")

    if not SKIP_PMT:
        print(f"  p-MT (N={N_pMT})...", end=" ", flush=True)
        m_pMT, s_pMT = measure(lfunc_pMT, (p, r, d, N_pMT, pre, dig), n_rep)
        print(f"{m_pMT:.5f} s")
    else:
        m_pMT, s_pMT = None, None

    print(f"  Diamond (N={N_Dia})...", end=" ", flush=True)
    m_Dia, s_Dia = measure(lfunc_Diamond, (p, d, N_Dia, pre, dig), n_rep)
    print(f"{m_Dia:.5f} s")

    rows[d] = (m_KL, s_KL, m_pMT, s_pMT, m_Dia, s_Dia)

print()
print(f"{'d':>4}  {'Alg. KL (s)':>14}  {'Alg. Diamond (s)':>16}", end="")
if not SKIP_PMT:
    print(f"  {'Alg. p-MT (s)':>14}", end="")
print()
print("-" * (4 + 2 + 14 + 2 + 16 + (2 + 14 if not SKIP_PMT else 0)))
for d in d_list:
    m_KL, s_KL, m_pMT, s_pMT, m_Dia, s_Dia = rows[d]
    line = f"{d:>4}  {m_KL:>14.5f}  {m_Dia:>16.5f}"
    if not SKIP_PMT:
        line += f"  {m_pMT:>14.5f}"
    print(line)
