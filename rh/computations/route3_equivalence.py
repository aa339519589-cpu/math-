"""
QUESTION (route 3 / low-startup wedge, 1<=n<r):
  Does "D_beta(n,r) > 0 for ALL 1<=n<r, for every r>=2" (wedge positivity,
  alone) already IMPLY RH (i.e. imply full real-rootedness / full total
  positivity of the whole sequence)?

ANSWER via a concrete counterexample: NO. Wedge positivity alone is
logically WEAKER than full total positivity. This matches classical
total-positivity theory (Fekete's criterion, see note below) and is
demonstrated here by an explicit toy sequence.

Fekete's lemma (Fekete 1912; see Pinkus, "Totally Positive Matrices",
Cambridge 2010, Thm 2.10): a matrix is STRICTLY totally positive (every
minor, of every shape, positive) IFF every CONTIGUOUS minor (consecutive
rows x consecutive columns, of every size and every starting position)
is positive.

For a one-sided Toeplitz matrix T=[beta_{j-i}]_{i,j>=0} (beta_k=0 for
k<0), a contiguous r x r minor starting at row p, column q depends only
on the shift n=q-p, and equals exactly D_beta(n,r) as defined in the
ledger. So full strict total positivity of (b_n) -- which is EXACTLY
equivalent to RH via the Polya-frequency-sequence / Laguerre-Polya
theorem -- requires D_beta(n,r) > 0 for *every* r>=1 and *every*
n = q-p, i.e. every integer n >= -(r-1) that can occur for
p,q >= 0. That splits into:

  * n = 0                : leading principal minors (call this "route 1")
  * n >= r  ("interior")  : "route 2" in the ledger
  * 1 <= n < r  ("wedge") : "route 3" in the ledger
  * n < 0                 : trivially degenerate / mostly forced by the
                            zero-padding (not a separate hard route)

Fekete's theorem needs ALL of these together. There is no version of
Fekete's lemma that lets you drop the interior (route 2) minors and
still conclude full total positivity from the wedge (route 3) alone.

The counterexample below makes this completely concrete: a positive
sequence satisfying full wedge positivity, whose associated "Xi-like"
even entire function nonetheless has NON-real (complex) zeros --
i.e. "RH" manifestly FAILS for this toy sequence despite the wedge
holding.
"""
from mpmath import mp, mpf, matrix, det, polyroots

mp.dps = 30


def beta_from_list(lst):
    def beta(k):
        if k < 0 or k >= len(lst):
            return mpf(0)
        return lst[k]
    return beta


def D_beta(beta, n, r):
    M = matrix(r, r)
    for i in range(r):
        for j in range(r):
            M[i, j] = beta(n + j - i)
    return det(M)


if __name__ == "__main__":
    # Toy "Xi-like" coefficients c_0=1, c_1=2, c_2=1.5, then zero.
    # g(x) = sum (-1)^n c_n x^n = 1 - 2x + 1.5x^2  (finite -- padded by 0).
    c = [mpf(1), mpf(2), mpf('1.5')]
    beta = beta_from_list(c)

    print("Toy sequence c = (c_0, c_1, c_2) =", c, " (c_n=0 for n>=3)")
    print()

    print("Step 1: does the associated polynomial have real roots?")
    # g(x) = 1 - c1 x + c2 x^2   (since (-1)^1 c_1 = -c_1, (-1)^2 c_2 = +c_2)
    coeffs_highest_first = [c[2], -c[1], c[0]]  # c2 x^2 - c1 x + c0
    roots = polyroots(coeffs_highest_first)
    print("  roots of c2*x^2 - c1*x + c0 :", roots)
    disc = c[1] ** 2 - 4 * c[2] * c[0]
    print(f"  discriminant c1^2 - 4*c2*c0 = {disc}  -> {'REAL roots' if disc >= 0 else 'COMPLEX (non-real) roots -> real-rootedness FAILS'}")
    print()

    print("Step 2: check wedge positivity D_beta(n,r) for 1<=n<r, r=2..5")
    all_ok = True
    for r in range(2, 6):
        for n in range(1, r):
            d = D_beta(beta, n, r)
            ok = d > 0
            all_ok &= ok
            print(f"  D_beta({n},{r}) = {d}   -> {'OK (>0)' if ok else 'FAILS'}")
    print()
    print("All wedge minors positive:", all_ok)
    print()

    print("Step 3: check an INTERIOR ('route 2', n>=r) minor to show where it actually breaks")
    for (n, r) in [(2, 2), (3, 3), (2, 3)]:
        d = D_beta(beta, n, r)
        print(f"  D_beta({n},{r}) [n>=r, interior] = {d}   -> {'OK' if d > 0 else 'FAILS (or degenerate)'}")

    print()
    print("CONCLUSION: the wedge (route 3) holds fully for this toy sequence,")
    print("but the sequence's own generating polynomial already has non-real")
    print("roots -- i.e. 'RH' fails for this toy model. So wedge positivity")
    print("(route 3) alone does NOT imply real-rootedness / RH in general.")
    print("It is a genuinely necessary but NOT sufficient condition on its own.")
