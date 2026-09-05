# Is Route 3 (low-startup wedge) alone equivalent to RH?

**Verdict: 不成立 — Route 3 alone is NOT logically equivalent to RH.**
This is proved below, not conjectured from numerics. The ledger's own
earlier caution ("这个结论不能在...尚未闭合时被认证", addendum §六) was
correct to withhold certification.

## I. 成立性裁决

### The structural reason (real theorem, not folklore)

RH ⟺ `(b_n)` is a Pólya-frequency sequence ⟺ the semi-infinite Toeplitz
matrix `T = [β_{j-i}]_{i,j≥0}` (β_k=0, k<0) is **totally positive**
(every minor, of every shape — not just contiguous ones — is ≥0).

**Fekete's lemma** (Fekete 1912; Pinkus, *Totally Positive Matrices*,
CUP 2010, Thm 2.10): a matrix is (strictly) totally positive **iff**
every *contiguous* minor (consecutive rows × consecutive columns, every
size, every starting position) is positive.

For a Toeplitz matrix, a contiguous `r×r` minor starting at row `p`,
column `q` depends only on the shift `n = q-p` and equals exactly
`D_β(n,r)`. So Fekete's criterion says: RH ⟺ `D_β(n,r) > 0` for
**every** `r ≥ 1` and **every achievable shift** `n`, which splits into
three disjoint families as `p,q` range over `{0,1,2,...}`:

| shift | ledger name | what it tests |
|---|---|---|
| `n = 0` | (route 1?) | leading principal minors |
| `n ≥ r` | route 2 | interior minors, no zero-padding involved |
| `1 ≤ n < r` | route 3 | wedge minors, touching the zero-padding edge |

Fekete's theorem needs **all three families simultaneously**. There is
no version of it that lets any one family stand in for the whole. So on
pure structural grounds, route 3 by itself can only ever be a
**necessary**, not sufficient, condition for RH.

## II. Explicit counterexample (not just structural argument)

To make "necessary but not sufficient" concrete rather than abstract,
here is an actual family of positive-coefficient sequences where the
**entire wedge holds over an arbitrarily long finite range, while
real-rootedness is exactly, provably false for every member of the
family**:

```
g_eps(x) = (1 + eps*x^2) * e^{-x},   eps > 0.
```

- Zeros of `g_eps`: `e^{-x}` is nowhere zero; `1+eps*x^2=0` at
  `x = ± i/sqrt(eps)` — **non-real for every eps>0, exactly, no
  approximation involved.** So the toy "RH" is false for every member
  of this family, unconditionally.
- Taylor coefficients `c_m = 1/m! + eps/(m-2)!` (second term absent for
  `m<2`) are manifestly **positive for every m** — the same baseline
  positivity that holds unconditionally for the real `b_n`.
- As `eps → 0+`, `g_eps → e^{-x}`, whose coefficients `1/n!` are a
  textbook Pólya-frequency (totally positive) sequence with **no**
  failures at any finite `(n,r)`. By continuity, for smaller `eps` the
  first wedge failure is pushed out to larger and larger `(n,r)`.

Computed (`computations/route3_equivalence4.py`, cross-checked at 50
and 200 significant digits to rule out numerical underflow artifacts):

| eps | first wedge failure `(n,r)` |
|---|---|
| 0.1   | `(2,4)` |
| 0.01  | `(5,9)` |
| 0.001 | `(15,21)`, value `≈ -8.2e-383` (genuine, not underflow — confirmed at 200 digits) |

So: **the wedge failure point can be pushed arbitrarily far out by
shrinking `eps`, while real-rootedness fails identically (and provably)
for every `eps>0`.** No finite computational check of the wedge —
however far extended — can ever certify real-rootedness, and there is
no generic mechanism by which "wedge holds for a long stretch" implies
"wedge holds forever," let alone "implies real-rootedness."

## III. Route decision

- **X (recommended):** keep route 3 as one-third of the Fekete
  decomposition. Proving `D_β(n,r)>0` for the true Ξ-sequence, for all
  `1≤n<r`, remains a legitimate, well-defined, worthwhile target — it
  is a genuine necessary condition (and, per the classical `r=2` Turán
  inequality, at least the smallest case is already unconditionally
  true for the real sequence). It should NOT be marketed, even if fully
  closed, as "= RH" on its own.
- **Y (only if route 1 and route 2 are already fully, independently
  closed):** if the project's route 1 (`n=0`) and route 2 (`n≥r`)
  results are correct and complete, then closing route 3 would
  complete the Fekete cover and **would** constitute a genuine proof of
  RH via this method. I have not independently re-audited routes 1/2's
  claimed proofs in this session — that re-audit is a prerequisite
  before route 3's closure could be promoted all the way to "RH proved."

## IV. Ledger update

- **已确认事实:** Fekete's lemma is the correct, citable reduction
  theorem connecting Toeplitz total positivity to contiguous minors;
  route 3 corresponds to exactly one of three necessary minor families.
- **成立点:** an explicit, provable (not numerical) counterexample
  family `g_eps` demonstrates wedge-positivity-over-an-arbitrarily-long-
  range does not imply real-rootedness in general.
- **失败 / 错误点:** any prior working assumption that "route 3 alone
  ⟺ RH" is now positively refuted in general, not merely unconfirmed.
- **仍属假设:** whether the *specific* Ξ-sequence's wedge holds for
  *all* r (not just the r≤6 numerically checked earlier) remains open;
  the counterexample says nothing about whether the real sequence
  happens to satisfy the full wedge — only that doing so would not, by
  itself, constitute a proof of RH.
- **未知 / 待验证:** the actual correctness of the claimed-closed route
  1 and route 2 results — not re-derived in this session, only assumed
  as given by the user.
- **下一步路线:** treat route 3 as a genuine but partial target; do not
  reopen it as an "RH ⟺" claim. If the team wants a real RH proof via
  this Fekete/total-positivity method, route 3 must be paired with a
  verified route 1 + route 2, not pursued as a standalone equivalence.
