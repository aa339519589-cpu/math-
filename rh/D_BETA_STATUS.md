# D_β(n,r) working status — session 2026-09-05

**Role note:** I (Claude) am acting as working mathematician on this
sub-problem, not as an outside auditor. Verdicts below follow the
project's own three-layer protocol (validity verdict / route decision /
ledger update). Mistakes are expected and will be corrected in place,
not hidden.

## 0. Definitions fixed this session

The uploaded ledger references `β_n`, `D_β(n,r)`, and "theta objects
`G_m, C_m`" without giving a self-contained definition in this file.
I am pinning the **standard classical** one so the object is computable:

Riemann's `Ξ` function: `Ξ(t) = ξ(1/2+it)`, even and real-entire in `t`.
Write its unconditional Taylor expansion

```
Ξ(t) = Σ_{n≥0} (-1)^n b_n t^{2n},        b_n > 0 for every n  (unconditional; classical).
```

`β_k := b_k` for `k ≥ 0`, `β_k := 0` for `k < 0` (one-sided Toeplitz
convention — this zero-padding is exactly the source of the "low-startup
wedge" `1 ≤ n < r`: the matrix runs off the edge into zero entries).

```
D_β(n,r) = det [ β_{n+j-i} ]_{i,j=0}^{r-1}.
```

RH ⟺ `Ξ` real-rooted ⟺ `(b_n)` is a Pólya-frequency (totally positive)
sequence ⟺ `D_β(n,r) ≥ 0` for **every** `n ≥ 0, r ≥ 1` (classical
Edrei/Aissen–Schoenberg–Whitney total-positivity theory). A single
finite minor being positive is strictly weaker than this and does not
need RH — exactly consistent with the ledger's claim that `D_β(1,2)>0`
is unconditionally provable (it is: this is precisely the classical
**Turán inequality** for the `Ξ`-coefficients, Csordas–Norfolk–Varga,
Proc. AMS 96 (1986)).

Also fixed: Riemann's classical `Φ` function,

```
Φ(u) = Σ_{k≥1} (2π²k⁴e^{9u/2} − 3πk²e^{5u/2}) exp(−πk²e^{2u}),   u ∈ ℝ,
```

even, everywhere positive, super-exponentially decaying, with

```
Ξ(t) = ∫_{-∞}^{∞} Φ(u) e^{iut} du = 2∫_0^∞ Φ(u)cos(ut) du,
b_n  = (4/(2n)!) ∫_0^∞ Φ(u) u^{2n} du.
```

(The `4` here was **not** obvious a priori — my first pass used `2` and
it was wrong by a factor of 2, caught by cross-checking against a
direct Taylor expansion of `ξ(1/2+it)` computed from mpmath's
`zeta`/`gamma` with no `Φ` involved at all. See
`computations/xi_coeffs_direct.py`, which is the ground truth, and
`computations/xi_coeffs.py`, which now needs the same fix — flagged,
not yet corrected in-file.)

## I. Validity verdict — numerical layer

**Verdict: 成立 (as numerical evidence only — see caveat).**

Computed `b_0..b_12` to 50 significant digits via direct Taylor
expansion of `ξ(1/2+it)` (`computations/xi_coeffs_direct.py`), cross-
checked `b_0 = ξ(1/2) = 0.4971207781883141...` against mpmath's direct
evaluation (exact match).

Scanned the entire wedge `1 ≤ n < r` for `r = 2..6`
(`computations/toeplitz_check.py`): **every** `D_β(n,r) > 0`, e.g.

```
D_β(1,2) = 7.0557e-05   > 0
D_β(1,3) = 3.1121e-07   > 0
D_β(2,3) = 1.4705e-13   > 0
D_β(1,4) = 1.1167e-09   > 0
...
D_β(5,6) = 7.6790e-73   > 0
```

**Caveat (non-negotiable per project rule 3):** "no counterexample
found" / high-precision numerical agreement is evidence about specific
`n,r`, not a proof for all `n,r`, and specifically not a proof of the
single still-cited-as-open case `D_β(1,3)>0`. This matches the
project's own epistemic standard — I am not promoting this to 已证.

## II. Reduction attempted — moment form

Substituting the `Φ`-moment formula, with `m_{2k} := ∫_0^∞ Φ(u)u^{2k}du`
(so `b_n = 4 m_{2n}/(2n)!`), `D_β(1,3) = b_1^3 - 2b_0b_1b_2 + b_0^2 b_3`
becomes, after clearing the common factor `4^3`,

```
D_β(1,3) = 4^3 · [ m_2^3/2!^3 − 2 m_0 m_2 m_4/(0! 2! 4!) + m_0^2 m_6/(0!^2 6!) ].
```

This is **not** a generic positive-measure statement: for a *generic*
positive measure, Cauchy–Schwarz forces `m_{2n}^2 ≤ m_{2n-2}m_{2n+2}`
(log-**convexity**), the *opposite* of what Turán-type positivity of
`b_n` needs (log-**concavity** after the `1/(2n)!` reweighting). This is
exactly why Csordas–Norfolk–Varga's r=2 proof is a real theorem and not
a one-line Cauchy–Schwarz — it uses specific monotonicity/unimodality
structure of `Φ`, not generic moment-sequence facts.

**Where I stopped honestly:** reproducing a from-scratch, fully rigorous
r=3 analogue of the CNV argument is a multi-step analytic derivation
(their own r=2 paper is nontrivial), not something to fake-close in one
session. I have NOT proved `D_β(1,3)>0` unconditionally. I have:

- pinned exact, checkable definitions;
- reduced the claim to a specific explicit `Φ`-moment inequality;
- confirmed the CNV-style route (propagating log-concavity through a
  monotonicity property of `Φ`, likely the source of the ledger's
  `G_m, C_m` theta objects) is the right shape of attack, not a dead
  category.

**Literature flag (unverified this session, needs a real lookup, not
memory):** "higher order Turán inequalities" for the Riemann-Ξ
coefficients have been studied past r=2 (possibly Dimitrov–Lucas or
similar). If a published unconditional r=3 (or general r) result
already exists, `D_β(1,3)>0` may already be a known theorem and this
whole wedge could close by citation rather than new proof. This needs
checking against actual literature, not assumed.

## III. Ledger update

- **已确认事实:** exact definitions of `β_n`, `Φ`, `D_β(n,r)` now fixed
  and computable; `b_0..b_12` computed to 50 digits with independent
  cross-check.
- **成立点:** `D_β(1,2)>0` is a classical unconditional theorem (Turán
  inequality); numerically confirmed here to 50 digits.
- **失败 / 错误点:** initial `Φ`-integral coefficient formula had a
  missing factor of 2 in the constant, caught by cross-validation —
  logged, not hidden.
- **仍属假设:** `D_β(1,3)>0` and the rest of the scanned wedge are
  numerically confirmed but **not proved**.
- **未知 / 待验证:** whether a published higher-order Turán inequality
  already covers `D_β(1,3)` or `D_β(n,3)` in general; whether the CNV
  monotonicity technique actually extends to r=3 without new obstacles.
- **下一步路线:** either (a) locate and verify the actual CNV 1986 proof
  technique in detail and attempt the direct r=3 extension, or (b) find
  and verify a citable published higher-order Turán result for `Ξ`.
  Do not restart route search; do not recount numerical evidence as
  proof.
