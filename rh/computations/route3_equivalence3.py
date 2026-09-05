"""
Third attempt, respecting the REAL constraint of the problem: c_n > 0 for
all n (this is a genuinely established, unconditional fact for the actual
Riemann Xi coefficients b_n -- not up for debate). The open question is
purely about the TOTAL-POSITIVITY / TOEPLITZ-MINOR sign pattern given
that baseline.

Take g(x) = (1+x^2) e^{-x}.
  - zeros of g: e^{-x} is nowhere zero; 1+x^2=0 at x=+-i -- PURELY
    IMAGINARY, i.e. non-real. So g is explicitly NOT real-rooted (the
    toy analogue of "RH fails").
  - Taylor coefficients: g(x) = sum (-1)^m c_m x^m with
        c_m = 1/m! + 1/(m-2)!   (second term absent for m<2),
    and c_m > 0 for EVERY m (both summands are nonnegative, matching the
    genuine unconditional b_n>0 fact for the real Xi function).

This is therefore a valid, apples-to-apples toy model: positive
coefficients, alternating-sign generating function, but manifestly
non-real-rooted. We now check whether the WEDGE (route 3, 1<=n<r) can
still hold across a wide range while the sequence is definitely not a
genuine Polya-frequency / totally positive sequence.
"""
from mpmath import mp, mpf, matrix, det, factorial

mp.dps = 40


def c(m):
    m = int(m)
    term1 = 1 / factorial(m)
    term2 = 1 / factorial(m - 2) if m >= 2 else mpf(0)
    return term1 + term2


def beta(k, n_terms=40):
    if k < 0:
        return mpf(0)
    return c(k)


def D_beta(n, r):
    M = matrix(r, r)
    for i in range(r):
        for j in range(r):
            M[i, j] = beta(n + j - i)
    return det(M)


if __name__ == "__main__":
    print("c_0..c_8 =", [float(c(m)) for m in range(9)])
    print("(all positive, as required)")
    print()
    print("g(x) = (1+x^2) e^{-x} has zeros at x = +-i  -- NOT real-rooted.")
    print("(toy 'RH' is FALSE for this sequence)")
    print()

    print("Wedge scan (route 3 shape), 1<=n<r, r=2..12:")
    wedge_fail = []
    for r in range(2, 13):
        for n in range(1, r):
            d = D_beta(n, r)
            if d <= 0:
                wedge_fail.append((n, r, d))
    if not wedge_fail:
        print("  ALL wedge minors D_beta(n,r), 1<=n<r, r=2..12, are POSITIVE.")
    else:
        print("  wedge fails at:", wedge_fail[:10])

    print()
    print("Interior scan (route 2 shape), n>=r, r=1..8, n=r..r+10:")
    interior_fail = []
    for r in range(1, 9):
        for n in range(r, r + 11):
            d = D_beta(n, r)
            if d <= 0:
                interior_fail.append((n, r, d))
    if not interior_fail:
        print("  ALL interior minors are POSITIVE too (would be a genuine surprise).")
    else:
        print(f"  interior FAILS at {len(interior_fail)} points, first few:")
        for f in interior_fail[:10]:
            print("   ", f)
