"""
Check the specific ledger claims:
   D_beta(1,2) > 0   (classical Turan inequality -- proved unconditionally,
                       Csordas-Norfolk-Varga 1986)
   D_beta(1,3) > 0   (the "low-startup wedge" case the ledger says was
                       still open / being supplemented)

Convention: beta_k = b_k for k >= 0, and beta_k = 0 for k < 0
(one-sided Toeplitz matrix built from a power series with only
nonnegative-index coefficients -- this IS the source of the "wedge"
difficulty: near n < r the matrix runs off the edge into zero entries).

D_beta(n,r) = det [ beta_{n+j-i} ]_{i,j=0}^{r-1}
"""
import sys
sys.path.insert(0, "/home/user/math-/rh/computations")
from xi_coeffs_direct import Xi
from mpmath import mp, mpf, taylor, matrix, det

mp.dps = 50

N_TERMS = 12
coeffs = taylor(lambda t: Xi(t).real, 0, 2 * N_TERMS)
b = [((-1) ** n) * coeffs[2 * n] for n in range(N_TERMS + 1)]


def beta(k):
    return b[k] if k >= 0 else mpf(0)


def D_beta(n, r):
    M = matrix(r, r)
    for i in range(r):
        for j in range(r):
            M[i, j] = beta(n + j - i)
    return det(M)


if __name__ == "__main__":
    print("b_n for n=0..%d:" % N_TERMS)
    for n, v in enumerate(b):
        print(f"  b_{n} = {v}")
    print()

    d12 = D_beta(1, 2)
    print(f"D_beta(1,2) = b_1^2 - b_0*b_2 = {d12}   -> {'>0 OK' if d12 > 0 else 'FAILS'}")

    d13 = D_beta(1, 3)
    print(f"D_beta(1,3) = {d13}   -> {'>0 OK' if d13 > 0 else 'FAILS'}")

    print()
    print("Wedge region 1<=n<r, scanning r=2..6, n=1..r-1:")
    for r in range(2, 7):
        for n in range(1, r):
            d = D_beta(n, r)
            print(f"  D_beta({n},{r}) = {d}   -> {'OK' if d > 0 else 'FAILS'}")
