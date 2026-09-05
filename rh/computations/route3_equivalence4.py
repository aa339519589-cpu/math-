"""
Fourth attempt: g_eps(x) = (1 + eps*x^2) * e^{-x}, eps>0 small.

  - zeros of g_eps: x = +- i/sqrt(eps)  -- ALWAYS non-real for any eps>0,
    however small. So "RH" (real-rootedness) fails for every eps>0,
    with no exception.

  - as eps -> 0+, g_eps -> e^{-x}, whose coefficient sequence 1/n! is a
    textbook Polya-frequency (totally positive) sequence. So by
    continuity, for eps small enough, any FIXED FINITE family of wedge
    minors D_beta(n,r) (r up to any fixed bound we choose to check)
    should stay positive, since they converge to the values for the
    (genuinely PF) 1/n! sequence as eps -> 0.

This is exactly the separation needed: for suitable small eps, wedge
positivity holds over a large *finite* range we can check computationally,
while real-rootedness/"RH" is *provably* false for that eps (roots
manifestly non-real) -- no computation needed for that part, it's exact.

Caveat, stated plainly: this shows wedge-positivity checked up to any
FIXED finite (n,r) range does not detect the failure -- it does not
(cannot, by finite computation) rule out that the FULL infinite wedge
(every r, no bound) might still somehow fail for this eps at some very
large r we didn't check. But it already refutes the naive claim
"wedge positive for a long stretch => real-rooted"; and it isolates
exactly why route 3 needs an honest infinite/asymptotic argument, not
finite-range numerics, to ever be promoted to an RH-equivalent theorem.
"""
from mpmath import mp, mpf, matrix, det, factorial

mp.dps = 50


def make_c(eps, n_terms):
    # g_eps(x) = (1+eps x^2) e^{-x} = e^{-x} + eps x^2 e^{-x}
    # [x^m] e^{-x} = (-1)^m/m!
    # [x^m] eps x^2 e^{-x} = eps * (-1)^{m-2}/(m-2)!  for m>=2
    c = []
    for m in range(n_terms + 1):
        t1 = mpf(1) / factorial(m)
        t2 = eps / factorial(m - 2) if m >= 2 else mpf(0)
        c.append(t1 + t2)  # matches (-1)^m sign convention automatically
    return c


def D_beta(c, n, r):
    def beta(k):
        return c[k] if 0 <= k < len(c) else mpf(0)
    M = matrix(r, r)
    for i in range(r):
        for j in range(r):
            M[i, j] = beta(n + j - i)
    return det(M)


if __name__ == "__main__":
    N_TERMS = 40
    for eps in [mpf('0.1'), mpf('0.01'), mpf('0.001')]:
        c = make_c(eps, N_TERMS)
        print(f"=== eps = {eps} ===  (zeros of g at x = +-i/sqrt(eps) = +-{1/eps**mpf('0.5')}i, non-real)")
        fails = []
        for r in range(2, 16):
            for n in range(1, r):
                d = D_beta(c, n, r)
                if d <= 0:
                    fails.append((n, r, float(d)))
        if not fails:
            print(f"  wedge D_beta(n,r), 1<=n<r, r=2..15: ALL POSITIVE")
            print(f"  yet real-rootedness is EXACTLY, PROVABLY false for this eps (roots +-i/sqrt(eps)).")
        else:
            print(f"  wedge FAILS at (first 5): {fails[:5]}")
        print()

    print("=== extended check at eps=0.001, r up to 25, plus interior minors ===")
    eps = mpf('0.001')
    c = make_c(eps, 60)
    wedge_fails = []
    for r in range(2, 26):
        for n in range(1, r):
            d = D_beta(c, n, r)
            if d <= 0:
                wedge_fails.append((n, r, float(d)))
    print("wedge fails (r up to 25):", wedge_fails if wedge_fails else "NONE -- all positive")

    interior_fails = []
    for r in range(1, 15):
        for n in range(r, r + 15):
            d = D_beta(c, n, r)
            if d <= 0:
                interior_fails.append((n, r, float(d)))
    print("interior fails (route-2 shape):", interior_fails[:10] if interior_fails else "NONE in range checked")
