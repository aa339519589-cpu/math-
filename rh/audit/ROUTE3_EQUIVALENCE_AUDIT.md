# Audit: the claimed proof of  R₃ ⟺ RH  (rectangular Schur route)

**Verdict: 条件成立.** The skeleton is sound and I could not break it.
It is **not** 成立, because three load-bearing steps are asserted rather
than proved. It is definitely **not** 不成立 — I found no error.

## 0. Retraction of my own previous session verdict

My earlier claimed disproof ("route 3 alone ⇏ RH", via
`g_ε(x)=(1+εx²)e^{−x}`) is **WITHDRAWN as out of scope.**

That family carries an essential `e^{−x}` factor, i.e. a genus-0 *order-1*
Hadamard factor. This argument's hypothesis class is
`A(z)=∏(1+x_ν z)`, `Σ|x_ν|<∞`, **no exponential factor** — and the real
Ξ genuinely lives there: with `z=−t²`, `Ξ` of order 1 in `t` becomes
**order 1/2 in z**, so no exponential factor is admissible and
`Σ|x_ν| = Σ|ρ−½|^{−2} < ∞` converges. My counterexample sat outside the
hypothesis class and refutes nothing here. Logged as an error, not
quietly dropped.

## I. What I verified as CORRECT

| # | claim | status |
|---|---|---|
| (1) | `D_β(n,r)/β_0^r = det[e_{n+j-i}] = s_{(r^n)}(X)` | ✅ correct — dual Jacobi–Trudi, conjugate of `(rⁿ)` is `(nʳ)` |
| — | genus-0 / `Σ|x_ν|<∞` hypothesis applies to the real Ξ | ✅ correct (order 1/2 in `z`) |
| §2 | RH ⟹ R₃ | ✅ correct (Schur SSYT monomial positivity) |
| (2) | `F_k(q)²−F_k(q−1)F_k(q+1)=F_{k−1}(q)F_{k+1}(q)` | ✅ verified numerically to 58 digits; genuine rectangular-Schur Q-system identity |
| — | exponentials cancel exactly when (2) is normalized by `(PRᵏ)^{2q}` | ✅ checked by hand |
| (4) | `d₀=d_m=0` | ✅ correct — the selection fills the rectangle completely, and `s_{(q^N)}(x₁..x_N)=(x₁⋯x_N)^q` exactly, so repeats produce no polynomial factor |
| (5) | `d₁=max_ν(ℓ_ν−1)` | ✅ confirmed in models (`h_q` of ℓ equal letters `= C(q+ℓ−1,ℓ−1)` ⟹ degree `ℓ−1`) |
| (6)→(7) | `d_{k−1}+d_{k+1}=2d_k−2`, `d₀=d_m=0` ⟹ `d_k=k(m−k)` | ✅ algebra correct |
| — | the `d_k=0` sub-case | ✅ clean: at top degree 0 the Haar step kills *every* `ω`, leaving `b_k=c+O(q^Cη^q)`, so the Turán difference decays exponentially while the RHS does not |
| §6 | `X₀={1,⅖e^{±iπ/3}}`: `D(1,r)>0 ∀r`, `D(2,3)=−776/15625` | ✅ **exact**: verified `D(2,3)=−0.049664`, and `D(1,r)>0` for `r≤44`; `p_m=1+2(⅖)^m cos(mπ/3)>0` and `exp` of a positive series ⟹ all `h_r>0` |

## II. Adversarial test — I tried to break it and failed

If the theorem is true, no letter multiset with positive real parts and a
non-real letter can satisfy the whole wedge. Searched ~70 configurations
(`audit/route3_falsify_search.py`): equal-letter blocks of size 1–4 with a
small non-real pair; a **doubled** non-real pair (the shape an actual
off-line Riemann zero quadruple produces); a non-real pair as top modulus;
near-real phases `θ=0.3, 0.1, 0.03`; a `1/k²` tail plus a tiny non-real
pair; plus **60 randomized configurations**.

**Result: every single one fails somewhere in the wedge. Zero survivors.**

Note also: unlike my withdrawn `e^{−x}` family, failures here appear at
*small* `(n,r)` — e.g. `X={1,εu,εū}` fails at `(n,r)=(2,3)` for every
`ε ∈ {0.4,0.2,0.1,0.05,0.02}` (value `∼ Cε³`). Removing the exponential
factor removes the "push the failure to infinity" phenomenon.

Caveat on the sharper reading: the *first* failure is not always at
`n=M+1` (e.g. the `1/k²` case first fails at `(11,12)` with `M+1=9`).
That is consistent — the mechanism's contradiction is asymptotic in `q`,
so it does not predict where the *smallest* failure sits.

## III. The three real gaps (why this is 条件成立, not 成立)

**G1 — equation (3) is asserted, not proved.** This is the technical
heart, and it gets one sentence ("再用 Cauchy–Binet，可得"). A real proof
needs: (a) an explicit spectral gap below `R` — true for Ξ, since the
moduli `|ρ−½|^{−2}` accumulate only at 0, so `{ν : |x_ν| ≥ R−δ}` is
finite, but this must be *stated and used*; (b) uniform control of the
infinite sub-dominant tail inside a fixed-size Hankel determinant;
(c) an actual proof that the coefficients are polynomials in `q`.

**G2 — the Haar-averaging descent is a multi-order induction compressed
into two sentences.** Worked out, the missing content is:
write `b_k(q)=c q^d + A(q)q^{d−1}+…` with `A` almost periodic. The cross
term in the Turán difference is `−c·q^{2d−1}·Δ²A(q)`, of degree
`2d−1 > 2d−2`, with oscillating sign (`Δ²A` has Haar mean 0), so
positivity forces `A` constant. One must then descend to degree `d−2`,
whose contribution lands at exactly `2d−2` and only perturbs the leading
coefficient to `c²d − cΔ²B(q)` — which cannot vanish identically, again
because `Δ²B` has mean 0. So **the claimed degree `2d_k−2` does survive**,
but via an induction the write-up does not contain. This is the gap
between "the conclusion is right" and "the proof is written".

**G3 — (5) needs an actual confluent-Vandermonde degree proof**, not
model verification.

**G4 — an unflagged arithmetic input.** `Re x_ν > 0` for every letter is
used at the end (to exclude `u=−1`). For ξ it holds (`|γ|≥14.13 > ½ ≥ |δ|`),
but it is arithmetic input, not formal consequence — it must be declared.

**G5 — non-effective.** Even fully proved, it bounds nothing about *where*
the wedge fails for the real Ξ.

## IV. Route decision — the strategic consequence

The question was: "要么闭环（路线3死了），要么不闭环（可以继续攻击）."

**Answer: 不闭环 — but you must NOT keep attacking route 3 as a cheaper
sub-target.** Those are different things, and conflating them is the
expensive mistake here.

- The equivalence argument survived every check I could run. Proceed on
  the working assumption that **R₃ ⟺ RH is true**.
- Therefore route 3 carries **the full difficulty of RH**. It is not the
  "leftover low-startup wedge after route 2". Any round spent on route 3
  in the belief that it is cheaper than RH is mispriced.
- What is genuinely open and genuinely smaller is **closing G1–G3**.
  That is a bounded, finite, checkable piece of work (asymptotics of
  fixed-size Hankel determinants with a spectral gap + an almost-periodic
  descent). It does not prove RH — it proves that route 3 is *not a
  shortcut*, which is exactly the information the project needs to stop
  re-investing in it.

**下一步路线 (recommended):** close G1 and G2 properly. If they close, the
ledger records `R₃ ⟺ RH` as 成立, and route 3 is retired as a separate
target — not because it is dead, but because it is RH itself wearing a
different notation.
