"""
KLpAdicL.sage — SageMath implementation of the Kubota-Leopoldt p-adic L-function.

Three algorithms are provided for computing L_p(r, chi_F) where chi_F is the
Kronecker character associated with a quadratic field F = Q(sqrt(d)):

  KL       — Kubota-Leopoldt original construction (partial sum of the
              defining limit formula for generalized Bernoulli numbers).
  p-MT     — p-adic Mellin transform (Riemann-sum approximation of the
              p-adic integral representation).
  Diamond  — Diamond's formula using Stirling numbers of the first kind;
              specialized to r=2 for improved speed.

The main entry point is the `KLpAdicL` class.  Module-level wrapper functions
are also provided for convenience.

Reference:
  Kimura, I. and Umeda, K., "Numerical Evaluation of the Kubota-Leopoldt
  p-adic L-Function at Positive Integers", JSIAM Letters (submitted).
"""

from functools import cache
from sage.misc.cachefunc import cached_method


class KLpAdicL:
    """Compute values of the Kubota-Leopoldt p-adic L-function.

    Parameters
    ----------
    p : int
        An odd prime.
    N : int
        Working p-adic precision (number of significant p-adic digits used
        internally).  Results are returned at this precision.

    Examples
    --------
    ::

        sage: kl = KLpAdicL(3, 20)
        sage: kl.lfunc_kl_original(2, 2, upto=1)
        2 + 3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^8 + 3^9 + ...
    """

    def __init__(self, p, N):
        """
        Parameters
        ----------
        p : int
            An odd prime.
        N : int
            Working p-adic precision.
        """
        self.p = p
        self.N = N
        self.Kp = pAdicField(p, N)

    @cached_method
    def chid_omegap_padic(self, d, r, h):
        """Return the Dirichlet character chi_d * omega_p^r as a p-adic character.

        Here chi_d is the Kronecker character (d/·) and omega_p is the
        Teichmüller character mod p (the unique character of order p-1).
        The generator of DirichletGroup(p) is raised to the h-th power to
        select the correct root of the (p-1)-th cyclotomic polynomial, so
        h must be coprime to p-1.

        Parameters
        ----------
        d : int
            A squarefree integer.  The conductor of chi_d is |d| or 4|d|
            according to whether d ≡ 1 (mod 4) or not.
        r : int
            Exponent of the Teichmüller character.
        h : int
            Must satisfy gcd(h, p-1) = 1.

        Returns
        -------
        DirichletCharacter
            The character chi_d * omega_p^r over pAdicField(p, N).
        """
        p = self.p
        if p == 2:
            raise NotImplementedError("p=2 is not supported.")
        if gcd(p - 1, h) != 1:
            raise ValueError("h must be coprime to p-1.")
        Kp = self.Kp
        DGp = DirichletGroup(p, Kp)
        if d % 4 == 2 or d % 4 == 3:
            e = abs(4 * d)
        else:
            e = abs(d)
        DGep = DirichletGroup(e * p, Kp)
        chid = kronecker_character(d)
        omegap = DGp.gen()^h
        chiomega = DGep(chid) * (DGep(omegap))^(r)
        return chiomega

    @staticmethod
    @cache
    def str1(n, m):
        """Return the (signed) Stirling number of the first kind s(n, m).

        Defined by the recurrence s(n, m) = s(n-1, m-1) - (n-1)*s(n-1, m)
        with s(n, n) = 1 and s(n, 0) = 0 for n >= 1.

        Examples
        --------
        ::

            sage: KLpAdicL.str1(2, 0)
            0
            sage: KLpAdicL.str1(3, 1)
            2
            sage: KLpAdicL.str1(5, 2)
            -50
        """
        if m == 0:
            return 0
        if n == m:
            return 1
        if m == 1:
            return (-1)^(n - 1) * factorial(n - 1)
        if m == n - 1:
            return -binomial(n, 2)
        if n == 0 and m == 0:
            return 0
        return KLpAdicL.str1(n - 1, m - 1) - (n - 1) * KLpAdicL.str1(n - 1, m)

    @staticmethod
    @cache
    def binomial2(n, x):
        """Return the n-th finite difference of 1/x, i.e. sum_{r=0}^n (-1)^r C(n,r)/(x+r).

        Examples
        --------
        ::

            sage: KLpAdicL.binomial2(0, x)
            1/x
            sage: KLpAdicL.binomial2(1, x)
            -1/(x + 1) + 1/x
        """
        return sum((binomial(n, r) * (-1)^r) / (x + r) for r in range(0, n + 1))

    @cached_method
    def lfunc_diamond(self, r, d, upto=None):
        """Compute L_p(r, chi_d) via Diamond's formula (Algorithm Diamond).

        Uses the series expansion involving Stirling numbers of the first kind
        and generalized binomial coefficients, valid for any positive integer r.

        Parameters
        ----------
        r : int
            A positive integer (the argument of the L-function).
        d : int
            A squarefree integer with gcd(d, p) = 1.
        upto : int, optional
            Truncation level N of the series.  Defaults to self.N + 1.

        Returns
        -------
        p-adic number
            Approximation of L_p(r, chi_d) at precision self.N.

        Examples
        --------
        ::

            sage: KLpAdicL(3, 20).lfunc_diamond(2, 2, upto=20)
            2 + 3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^8 + 3^9 + ...
        """
        p = self.p
        N = self.N
        if upto is None:
            upto = self.N + 1
        if is_squarefree(d) == False:
            return "Error:d must be squarefree"
        if d % 4 == 2 or d % 4 == 3:
            d = 4 * d
        chiomega = self.chid_omegap_padic(d, r - 1, 1)
        f_tilde = chiomega.conductor()
        total = 0
        for a in range(1, (p * f_tilde) + 1):
            if gcd(a, p) != 1:
                continue
            inner = 0
            for n in range(r - 2, upto):
                inner += ((-1)^(n + r)) / (factorial(n + 1)) * KLpAdicL.str1(n + 1, r - 1) * KLpAdicL.binomial2(n, a / (p * f_tilde))
            total += chiomega(a) * inner
        return (p * f_tilde)^(-r) * total

    def lfunc_diamond_2_chi(self, d, upto=None):
        """Compute L_p(2, chi_d) via Diamond's formula specialized to r=2 (Algorithm Diamond, r=2).

        This is an optimized form of lfunc_diamond for r=2, using factorial-based
        series coefficients that avoid computing Stirling numbers.

        Parameters
        ----------
        d : int
            A squarefree integer with gcd(d, p) = 1.
        upto : int, optional
            Truncation level N.  Defaults to self.N + 1.

        Returns
        -------
        p-adic number
            Approximation of L_p(2, chi_d) at precision self.N.

        Examples
        --------
        ::

            sage: KLpAdicL(3, 20).lfunc_diamond_2_chi(2, upto=20)
            2 + 3 + 3^2 + 3^3 + 2*3^4 + 2*3^5 + 2*3^6 + 2*3^8 + 3^9 + ...
        """
        if upto is None:
            upto = self.N + 1
        p = self.p
        N = self.N
        if not is_squarefree(d):
            return "Error:d must be squarefree"
        if d % 4 == 2 or d % 4 == 3:
            d = 4 * d
        chiomega = self.chid_omegap_padic(d, 1, 1)
        f_tilde = chiomega.conductor()
        return (sum(chiomega(a)*sum(factorial(n)*(p*f_tilde)^(n-1)/((n+1)*(product(a+i*p*f_tilde for i in range(n+1)))) for n in range(upto)) for a in range(1, (p*(f_tilde))+1) if (gcd(a,p)==1)))

    @staticmethod
    def _extract_pd_from_chi(chi):
        if isinstance(chi, (tuple, list)) and len(chi) == 2:
            return chi[0], chi[1]
        if isinstance(chi, dict) and "p" in chi and "d" in chi:
            return chi["p"], chi["d"]
        if hasattr(chi, "p") and hasattr(chi, "d"):
            return chi.p, chi.d
        raise ValueError("chi must provide p and d (e.g. (p, d) or {'p': p, 'd': d}).")

    def kl_padic_l(self, s, chi, upto=None, algorithm='default'):
        """Compute a value of the Kubota-Leopoldt p-adic L-function.

        Unified entry point that dispatches to the three algorithms.

        Parameters
        ----------
        s : int or p-adic number
            The argument of the L-function.
        chi : tuple, dict, or object
            Specification of the character.  Accepted forms:
              - (p, d)          — a 2-tuple
              - {'p': p, 'd': d} — a dict
              - any object with .p and .d attributes
        upto : int, optional
            Truncation level N.  Defaults to self.N + 1.
        algorithm : str, optional
            One of 'default', 'padicint', 'original', 'diamond'.
            'default' uses Algorithm Diamond for positive integer s and
            Algorithm p-MT otherwise.

        Returns
        -------
        p-adic number
            Approximation of L_p(s, chi) at precision self.N.
        """
        if upto is None:
            upto = self.N + 1
        _, d = self._extract_pd_from_chi(chi)
        if algorithm == 'padicint':
            return self.lfunc_padic_integral(s, d, upto)
        if algorithm == 'original':
            return self.lfunc_kl_original(s, d, upto)
        # algorithm == 'default'
        if s == 2:
            return self.lfunc_diamond_2_chi(d, upto)
        if (s in ZZ) and (s > 0):
            return self.lfunc_diamond(s, d, upto)
        return self.lfunc_padic_integral(s, d, upto)

    @staticmethod
    @cache
    def coprime_element(n):
        """Return the smallest integer >= 2 that is coprime to n."""
        a = 2
        while gcd(a, n) != 1:
            a += 1
        return a

    def lfunc_padic_integral(self, s, d, upto=None):
        """Compute L_p(s, chi_d) via the p-adic integral (Bernoulli measure) method (Algorithm p-MT).

        Approximates the p-adic Mellin transform integral by a finite Riemann sum.
        Works for any integer s (including non-positive integers).

        Parameters
        ----------
        s : int
            The argument of the L-function.
        d : int
            A squarefree integer with gcd(d, p) = 1.
        upto : int, optional
            Truncation level N (number of p-adic digits of integration mesh).
            Defaults to int(self.N / 2) + 1.

        Returns
        -------
        p-adic number
            Approximation of L_p(s, chi_d) at precision self.N.

        Examples
        --------
        ::

            sage: KLpAdicL(3, 12).lfunc_padic_integral(0, 1) == 3^-1 + 3^10 + O(3^18)
            True
            sage: KLpAdicL(3, 12).lfunc_padic_integral(-1, 1) == 2*3^-1 + 1 + 3 + ... + O(3^18)
            True
        """
        p = self.p
        if upto is None:
            upto = int(self.N / 2) + 1
        chiomega = self.chid_omegap_padic(d, -1, 1)
        chi = kronecker_character(d)
        K = self.Kp

        pN  = p^upto
        dpN = d * pN
        dp  = d * p
        alpha = KLpAdicL.coprime_element(dpN)

        S = IntegerModRing(dpN)
        inv_alpha = ZZ(S(alpha)^(-1))

        # Compute powers of the Teichmüller lift outside the sum to avoid
        # repeated p-adic exponentiations.
        tei_alpha       = K.teichmuller(alpha)^(1 - s)
        const_coeff     = -tei_alpha / (tei_alpha - chi(alpha) * K(alpha)^(1 - s))
        const_half      = K((alpha - 1) / 2)
        const_alpha_dpN = K(alpha) / K(dpN)

        total = sum(
            chiomega(a)
            * (K.teichmuller(a + pN) / K(a + pN))^s
            * (const_half + const_alpha_dpN * (a * inv_alpha - S(a * inv_alpha).lift()))
            for a in range(dpN) if gcd(a, dp) == 1
        )
        return const_coeff * total

    @cached_method
    def lfunc_kl_original(self, r, d, upto=None):
        """Compute L_p(r, chi_d) via the Kubota-Leopoldt original construction (Algorithm KL).

        Uses the limit formula derived from Witt's theorem for generalized
        Bernoulli numbers:

            L_p(r, chi) = (1/(r-1)) * lim_{N->inf} (1/(f*p^N))
                          * sum_{a=1, gcd(a,p)=1}^{f*p^N} chi*omega^{r-1}(a) * a^{1-r}

        Parameters
        ----------
        r : int
            A positive integer (the argument of the L-function).
        d : int
            A squarefree integer with gcd(d, p) = 1.
        upto : int, optional
            Truncation level N.  Defaults to int(self.N / 2) + 1.

        Returns
        -------
        p-adic number
            Approximation of L_p(r, chi_d) at precision self.N.

        Examples
        --------
        ::

            sage: KLpAdicL(3, 10).lfunc_kl_original(0, 1) == 3^-1 + O(3^9)
            True
            sage: KLpAdicL(3, 15).lfunc_kl_original(-1, 1) - 1/6 == O(3^4)
            True
        """
        if upto is None:
            upto = int(self.N / 2) + 1
        p = self.p
        N = upto
        chiomega = self.chid_omegap_padic(d, r - 1, 1)
        K = self.Kp

        pN  = p^upto
        dpN = d * p * pN   # = d * p^(upto+1)
        exp = 1 - r

        # Compute a^exp as an integer first when exp >= 0, then cast to K,
        # reducing the number of p-adic multiplications.
        if exp >= 0:
            total = sum(chiomega(a) * K(a^exp)
                        for a in range(1, dpN + 1) if a % p != 0)
        else:
            total = sum(chiomega(a) * K(a)^exp
                        for a in range(1, dpN + 1) if a % p != 0)

        return total / K((r - 1) * dpN)

    def kl_padic_l_diamond(self, s, chi):
        """Compute L_p(s, chi) using Algorithm Diamond.

        Thin wrapper around lfunc_diamond that accepts a chi specification
        in the same form as kl_padic_l.
        """
        _, d = self._extract_pd_from_chi(chi)
        return self.lfunc_diamond(s, d)


# ---------------------------------------------------------------------------
# Module-level convenience wrappers
# ---------------------------------------------------------------------------

_KL_INSTANCE = KLpAdicL(3, 20)


@cache
def chid_omegap_padic(d, p, r, h, N):
    """Module-level wrapper: return chi_d * omega_p^r over pAdicField(p, N)."""
    return KLpAdicL(p, N).chid_omegap_padic(d, r, h)


@cache
def str1(n, m):
    """Module-level wrapper: Stirling number of the first kind s(n, m)."""
    return KLpAdicL.str1(n, m)


@cache
def binomial2(n, x):
    """Module-level wrapper: sum_{r=0}^n (-1)^r C(n,r)/(x+r)."""
    return KLpAdicL.binomial2(n, x)


@cache
def lfunc_diamond(p, r, d, N):
    """Module-level wrapper: Algorithm Diamond for L_p(r, chi_d)."""
    return KLpAdicL(p, N).lfunc_diamond(r, d)


@cache
def lfunc_diamond_2_chi(p, d, N):
    """Module-level wrapper: Algorithm Diamond (r=2) for L_p(2, chi_d)."""
    return KLpAdicL(p, N).lfunc_diamond_2_chi(d)


@cache
def lfunc_kl_original(p, r, d, N):
    """Module-level wrapper: Algorithm KL for L_p(r, chi_d)."""
    return KLpAdicL(p, N).lfunc_kl_original(r, d)


def coprime_element(n):
    """Return the smallest integer >= 2 coprime to n."""
    return KLpAdicL.coprime_element(n)


def lfunc_padic_integral(p, s, d, N):
    """Module-level wrapper: Algorithm p-MT for L_p(s, chi_d)."""
    return KLpAdicL(p, N).lfunc_padic_integral(s, d)
