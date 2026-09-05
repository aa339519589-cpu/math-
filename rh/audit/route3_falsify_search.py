"""
Adversarial search: try to BREAK the claimed theorem  R_3 => RH.

If the theorem is true, then ANY letter multiset X with
  - Sum|x| < infinity (finite here),
  - all letters having positive real part,
  - at least one NON-REAL letter,
must have some wedge minor D(n,r) <= 0 with 1 <= n < r.

We search many configurations and record where the first failure occurs.
A configuration surviving the whole scanned wedge would falsify the theorem
(or push it beyond effective reach).

We also check the mechanism's sharper prediction: the failure should show up
in the family n = M+1, where M = #{letters of modulus > R}, R = max modulus
among non-real letters.
"""
from mpmath import mp, mpf, mpc, matrix, det, exp, pi, j as I
import random

mp.dps = 50
random.seed(20260905)


def elementary(letters, kmax):
    e = [mpc(0)] * (kmax + 1)
    e[0] = mpc(1)
    for x in letters:
        for k in range(min(kmax, len(letters)), 0, -1):
            e[k] = e[k] + e[k - 1] * x
    return e


def D(e, n, r):
    def E(k):
        return e[k] if 0 <= k < len(e) else mpc(0)
    M = matrix(r, r)
    for i in range(r):
        for j in range(r):
            M[i, j] = E(n + j - i)
    return det(M)


def analyze(X, label, rmax=26):
    e = elementary(X, 2 * rmax + 6)
    nonreal = [x for x in X if abs(x.imag) > mpf(10) ** -30]
    R = max(abs(x) for x in nonreal)
    M = sum(1 for x in X if abs(x) > R + mpf(10) ** -30)
    first = None
    for r in range(2, rmax + 1):
        for n in range(1, r):
            v = D(e, n, r).real
            if v <= 0:
                first = (n, r, v)
                break
        if first:
            break
    pred = M + 1
    if first:
        tag = "matches n=M+1" if first[0] == pred else f"n={first[0]} vs predicted M+1={pred}"
        print(f"  {label:<44} FAILS at (n,r)=({first[0]},{first[1]})  [{tag}]")
    else:
        print(f"  {label:<44} *** NO WEDGE FAILURE up to r={rmax} -- POSSIBLE COUNTEREXAMPLE ***")
    return first


print("=" * 78)
print("Adversarial search for a counterexample to  R_3 => RH")
print("=" * 78)

u60 = exp(I * pi / 3)

configs = []
# non-real pair well below a block of equal real letters
for N in [1, 2, 3, 4]:
    configs.append(([mpc(1)] * N + [mpf('0.1') * u60, mpf('0.1') * u60.conjugate()],
                    f"{N} x (1.0)  +  0.1*e^{{+-i pi/3}}"))
# repeated non-real pair (mimics an off-line Riemann zero quadruple)
configs.append(([mpc(1), mpf('0.3') * u60, mpf('0.3') * u60, mpf('0.3') * u60.conjugate(),
                 mpf('0.3') * u60.conjugate()], "1.0 + DOUBLED 0.3*e^{+-i pi/3}"))
# non-real pair as the TOP modulus
configs.append(([mpf('0.9') * u60, mpf('0.9') * u60.conjugate(), mpf('0.4'), mpf('0.15')],
                "non-real pair is TOP modulus"))
# tiny phase (nearly real non-real pair)
for th in ['0.3', '0.1', '0.03']:
    ph = exp(I * mpf(th))
    configs.append(([mpc(1), mpf('0.5') * ph, mpf('0.5') * ph.conjugate()],
                    f"1.0 + 0.5*e^{{+-i*{th}}}  (nearly real)"))
# decaying real tail plus a small non-real pair (closest to the Xi shape)
tail = [mpf(1) / (k * k) for k in range(1, 9)]
ph = exp(I * mpf('0.4'))
configs.append((tail + [mpf('0.013') * ph, mpf('0.013') * ph.conjugate()],
                "1/k^2 tail (k<=8) + tiny non-real pair"))

for X, label in configs:
    analyze(X, label)

print()
print("Randomized sweep (60 random configurations):")
survivors = 0
for trial in range(60):
    nreal = random.randint(0, 4)
    X = [mpf(random.uniform(0.2, 2.0)) for _ in range(nreal)]
    th = mpf(random.uniform(0.05, 3.0))
    rad = mpf(random.uniform(0.05, 1.5))
    ph = exp(I * th)
    X += [rad * ph, rad * ph.conjugate()]
    if all(x.real > 0 for x in X):
        e = elementary(X, 60)
        found = False
        for r in range(2, 25):
            for n in range(1, r):
                if D(e, n, r).real <= 0:
                    found = True
                    break
            if found:
                break
        if not found:
            survivors += 1
            print(f"  *** SURVIVOR (trial {trial}): {[mp.nstr(x, 5) for x in X]}")
print(f"  survivors out of 60: {survivors}")
