"""
Retry with a genuinely INFINITE (non-terminating) sequence -- the first
attempt used a finite (eventually-zero) sequence, which degenerates the
wedge itself at large r for a different, uninteresting reason (running
out of nonzero terms), so it did not actually test the real question.

Take g(x) = (1 + x + x^2) / (1 - r x),  0 < r < 1.

Numerator 1+x+x^2 has genuinely complex roots (primitive cube roots of
unity, scaled) -- so g, as an entire-function candidate for the ASW/
Laguerre-Polya class, is built from a "bad" (non-real-rooted) factor on
top of a "good" (single real pole) factor. Generic expectation: the
resulting coefficient sequence should FAIL full total positivity
(fail somewhere, wedge or interior), because one non-real-rooted factor
generically cannot be absorbed/cancelled by a single real pole.

We track SIGN of (-1)^n * [x^n] g(x) matching the beta_n>0 convention:
define c_n via g(x) = sum (-1)^n c_n x^n, i.e. c_n = (-1)^n [x^n] g(x).
"""
from mpmath import mp, mpf, matrix, det, taylor

mp.dps = 40


def make_beta(r_pole, n_terms=40):
    # g(x) = (1+x+x^2)/(1-r_pole x); Taylor coefficients via direct recursion:
    # (1 - r x) * g(x) = 1 + x + x^2  =>  g_n - r*g_{n-1} = [x^n](1+x+x^2)
    num = [mpf(1), mpf(1), mpf(1)]  # coefficients of 1+x+x^2
    g = [mpf(0)] * (n_terms + 1)
    for n in range(n_terms + 1):
        rhs = num[n] if n < len(num) else mpf(0)
        prev = g[n - 1] if n >= 1 else mpf(0)
        g[n] = rhs + r_pole * prev
    c = [((-1) ** n) * g[n] for n in range(n_terms + 1)]

    def beta(k):
        if k < 0 or k >= len(c):
            return mpf(0)
        return c[k]
    return beta, c


def D_beta(beta, n, r):
    M = matrix(r, r)
    for i in range(r):
        for j in range(r):
            M[i, j] = beta(n + j - i)
    return det(M)


if __name__ == "__main__":
    for r_pole in [mpf('0.3'), mpf('0.5'), mpf('0.7')]:
        beta, c = make_beta(r_pole, n_terms=30)
        print(f"=== pole r={r_pole}: c_0..c_6 = {[float(x) for x in c[:7]]} ===")
        wedge_ok = True
        first_wedge_fail = None
        for rr in range(2, 12):
            for n in range(1, rr):
                d = D_beta(beta, n, rr)
                if d <= 0:
                    wedge_ok = False
                    if first_wedge_fail is None:
                        first_wedge_fail = (n, rr, d)
        interior_ok = True
        first_interior_fail = None
        for rr in range(1, 8):
            for n in range(rr, rr + 15):
                d = D_beta(beta, n, rr)
                if d <= 0:
                    interior_ok = False
                    if first_interior_fail is None:
                        first_interior_fail = (n, rr, d)
        print(f"  wedge (1<=n<r, r=2..11) all positive? {wedge_ok}"
              + ("" if wedge_ok else f"  first fail at {first_wedge_fail}"))
        print(f"  interior (n>=r, r=1..7, n..n+15) all positive? {interior_ok}"
              + ("" if interior_ok else f"  first fail at {first_interior_fail}"))
        print()
