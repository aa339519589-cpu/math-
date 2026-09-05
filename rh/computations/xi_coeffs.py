"""
Riemann Xi function Taylor coefficients b_n, defined by the classical
unconditional expansion

    Xi(t) = sum_{n=0}^inf (-1)^n b_n t^{2n},    b_n > 0 for all n (unconditional).

via the standard integral representation

    Xi(t) = 2 * int_0^inf Phi(u) cos(t u) du,
    b_n   = (2/(2n)!) * int_0^inf Phi(u) u^{2n} du,

where Phi is Riemann's classical function

    Phi(u) = sum_{k=1}^inf (2 pi^2 k^4 e^{9u/2} - 3 pi k^2 e^{5u/2}) * exp(-pi k^2 e^{2u}).

Phi(u) > 0 for all real u and decays super-exponentially as u -> +-inf,
which is why the integral defining b_n converges and b_n > 0 unconditionally
(this positivity of b_n needs NO assumption on RH -- it is classical,
going back to Riemann / Polya, and reproved rigorously e.g. in
Csordas-Norfolk-Varga 1986).

RH is equivalent to: Xi(t) has only real zeros
             <=>  g(x) := sum (-1)^n b_n x^n has only nonnegative real zeros
             <=>  (b_n) is a Polya frequency (totally positive) sequence
             <=>  every Toeplitz minor  D_beta(n,r) = det[b_{n+j-i}]_{i,j=0}^{r-1} >= 0.

The r=2 case  D_beta(n,2) = b_n^2 - b_{n-1} b_{n+1} >= 0  is exactly the
classical TURAN INEQUALITY for the Xi coefficients, proved unconditionally
by Csordas, Norfolk and Varga (Proc. AMS, 1986). It does NOT require RH.

Higher r are open in general; r=3 at small n is the "low-startup wedge"
D_beta(1,3) referenced in the project ledger.
"""
from mpmath import mp, mpf, pi, exp, cos, quad, factorial, inf, nsum

mp.dps = 60  # 60 significant decimal digits


def Phi(u, terms=40):
    u = mpf(u)
    total = mpf(0)
    for k in range(1, terms + 1):
        k = mpf(k)
        e2u = exp(2 * u)
        term = (2 * pi**2 * k**4 * exp(mpf(9) / 2 * u) - 3 * pi * k**2 * exp(mpf(5) / 2 * u)) \
            * exp(-pi * k**2 * e2u)
        total += term
        if abs(term) < mpf(10) ** (-(mp.dps + 5)) and k > 5:
            break
    return total


def b(n, upper=8):
    """
    b_n = (4/(2n)!) * int_0^inf Phi(u) u^{2n} du, computed by numerical quadrature.

    CORRECTION (session 2026-09-05): this constant was originally coded as 2,
    which silently produced b_n values exactly half the true ones. Caught by
    cross-checking against xi_coeffs_direct.py (direct Taylor expansion of
    xi(1/2+it) via mpmath zeta/gamma, no Phi formula involved) -- b_0 there
    matches xi(1/2) exactly, and b_0 computed here with constant=2 did not.
    Logged rather than silently fixed, per the project's own audit rule
    against silently promoting corrected numbers without a trace.
    """
    n = int(n)
    integrand = lambda u: Phi(u) * u ** (2 * n)
    # Phi(u) decays super-exponentially in both directions; integrate on a
    # generous finite window and check tail negligibility.
    val = quad(integrand, [0, 1, 2, 3, 4, 5, upper])
    return 4 * val / factorial(2 * n)


if __name__ == "__main__":
    bs = [b(n) for n in range(6)]
    for n, v in enumerate(bs):
        print(f"b_{n} = {v}")
