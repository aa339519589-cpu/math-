"""
Cross-check: compute Riemann Xi(t) = xi(1/2+it) DIRECTLY from mpmath's
zeta/gamma, with no intermediate Phi-function formula, then read off the
Taylor coefficients b_n via  Xi(t) = sum (-1)^n b_n t^{2n}.

This is the ground truth used to validate (or catch bugs in) the
Phi-integral route in xi_coeffs.py.
"""
from mpmath import mp, mpf, mpc, pi, gamma, zeta, taylor, chop

mp.dps = 50


def xi(s):
    s = mpc(s)
    return 0.5 * s * (s - 1) * pi ** (-s / 2) * gamma(s / 2) * zeta(s)


def Xi(t):
    return xi(mpf(1) / 2 + 1j * mpf(t))


if __name__ == "__main__":
    # Taylor-expand Xi(t) about t=0 directly (Xi is real-analytic and even in t).
    coeffs = taylor(lambda t: Xi(t).real, 0, 10)
    print("Direct Taylor coefficients of Xi(t) about t=0 (should vanish at odd order):")
    for n, c in enumerate(coeffs):
        print(f"  t^{n}: {c}")

    print()
    print("b_n := (-1)^n * [t^{2n}] Xi(t):")
    for n in range(6):
        c2n = coeffs[2 * n]
        b_n = ((-1) ** n) * c2n
        print(f"  b_{n} = {b_n}")

    print()
    print("Cross-check b_0 against xi(1/2):", xi(mpf(1) / 2).real)
