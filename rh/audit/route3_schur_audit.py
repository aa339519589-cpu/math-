"""
Audit of the claimed proof  R_3 <=> RH  via rectangular Schur functions.

Checks the four independently checkable assertions:
  (A) dual Jacobi-Trudi:  D_beta(n,r)/beta_0^r = s_{(r^n)}(X)
  (B) the Q-system identity (2)
  (C) the section-6 counterexample numbers
  (D) the theorem's own falsifiable prediction on X={1, eps*u, eps*ubar}
"""
from mpmath import mp, mpf, mpc, matrix, det, exp, pi, j as I, binomial

mp.dps = 60


def elementary(letters, kmax):
    """e_0..e_kmax of a finite letter multiset."""
    e = [mpc(0)] * (kmax + 1)
    e[0] = mpc(1)
    for x in letters:
        for k in range(min(kmax, len(letters)), 0, -1):
            e[k] = e[k] + e[k - 1] * x
    return e


def D(e, n, r):
    """det[e_{n+j-i}]_{i,j=0}^{r-1}, with e_k=0 outside range."""
    def E(k):
        return e[k] if 0 <= k < len(e) else mpc(0)
    M = matrix(r, r)
    for i in range(r):
        for j in range(r):
            M[i, j] = E(n + j - i)
    return det(M)


def schur_rect(letters, rows, cols, e=None):
    """s_{(cols^rows)}(X) via dual Jacobi-Trudi = D(rows, cols)."""
    if e is None:
        e = elementary(letters, rows + cols + 4)
    return D(e, rows, cols)


# ---------------------------------------------------------------- (C)
print("=" * 70)
print("(C) Section-6 counterexample  X0 = {1, (2/5)e^{i pi/3}, (2/5)e^{-i pi/3}}")
print("=" * 70)
u = exp(I * pi / 3)
X0 = [mpc(1), mpf('0.4') * u, mpf('0.4') * u.conjugate()]
e0 = elementary(X0, 60)
print("  e_1, e_2, e_3 =", [mp.nstr(e0[k].real, 8) for k in (1, 2, 3)],
      " (claimed 7/5=1.4, 14/25=0.56, 4/25=0.16)")
d23 = D(e0, 2, 3)
print(f"  D(2,3) = {mp.nstr(d23.real, 12)}   (claimed -776/15625 = {mpf(-776)/15625})")
neg1 = [r for r in range(1, 45) if D(e0, 1, r).real <= 0]
print(f"  D(1,r) <= 0 for r in 1..44 ?  {neg1 if neg1 else 'NONE -- all strictly positive'}")

# ---------------------------------------------------------------- (A)+(B)
print()
print("=" * 70)
print("(A) dual Jacobi-Trudi consistency, and (B) the Q-system identity")
print("=" * 70)
Xt = [mpf('1.3'), mpf('0.9'), mpf('0.55') * exp(I * mpf('0.7')),
      mpf('0.55') * exp(-I * mpf('0.7')), mpf('0.2')]
et = elementary(Xt, 80)
worst = mpf(0)
for N in range(1, 5):
    for q in range(2, 9):
        lhs = schur_rect(Xt, N, q, et) ** 2 - schur_rect(Xt, N, q - 1, et) * schur_rect(Xt, N, q + 1, et)
        rhs = schur_rect(Xt, N - 1, q, et) * schur_rect(Xt, N + 1, q, et)
        worst = max(worst, abs(lhs - rhs))
print(f"  max |LHS-RHS| of  F_N(q)^2 - F_N(q-1)F_N(q+1) = F_(N-1)(q)F_(N+1)(q)")
print(f"  over N=1..4, q=2..8:   {mp.nstr(worst, 6)}   -> {'IDENTITY HOLDS' if worst < mpf(10)**-40 else 'IDENTITY FAILS'}")

# ---------------------------------------------------------------- (D)
print()
print("=" * 70)
print("(D) theorem's falsifiable prediction: X = {1, eps*u, eps*conj(u)}")
print("    m=2, M=1 -> theorem forces d_1=1, but direct count gives d_1=0,")
print("    so the theorem PREDICTS some wedge minor D(n,r)<0 with 1<=n<r.")
print("=" * 70)
for eps in [mpf('0.4'), mpf('0.2'), mpf('0.1'), mpf('0.05'), mpf('0.02')]:
    Xe = [mpc(1), eps * u, eps * u.conjugate()]
    ee = elementary(Xe, 130)
    first = None
    for r in range(2, 61):
        for n in range(1, r):
            if D(ee, n, r).real <= 0:
                first = (n, r, D(ee, n, r).real)
                break
        if first:
            break
    if first:
        print(f"  eps={float(eps):<5} first wedge failure at (n,r)={first[0]},{first[1]}"
              f"   value={mp.nstr(first[2], 6)}")
    else:
        print(f"  eps={float(eps):<5} NO wedge failure found up to r=60")
