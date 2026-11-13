# Sign structure, diagonal similarity, and the off-diagonal norm `τ`

**Context.** The certificate in `CalcChemistry/ComplementBound.lean` splits the complement
Hamiltonian as `H_C = D + O` with `D` diagonal, `D_i ≥ d`, and consumes a single scalar
`τ` with

> `|Q_O(x)| ≤ τ ‖x‖²` for all `x`,  where `Q_O(x) = Σ_i Σ_j O_ij x_i x_j`.

The question addressed here: *can sign structure — `O_ij` of mixed sign — certify a smaller
`τ` than the magnitude-only (Gershgorin/Schur) bound, e.g. by minimising `‖S O S⁻¹‖₂` over
positive diagonal `S`?*

**Answers, up front.**

1. For **symmetric** `O` the 2-norm route is a dead end: `inf_S ‖S O S⁻¹‖₂ = ‖O‖₂`
   exactly (Proposition 2.1). Similarity preserves the spectrum, the spectral norm
   dominates the spectral radius, and `‖O‖₂ = ρ(O)` for symmetric `O`. Moreover the
   constant `‖O‖₂` is *attained* by an eigenvector, so **no** sign argument can certify a
   smaller `τ` inside this formulation. The classical "diagonal scaling collapses the norm"
   phenomenon is real but is a statement about **non-normal** matrices
   (Bauer–Fike, Varah, Osborne) and about the **∞-norm**, not the 2-norm of a symmetric
   matrix.
2. The exact optimisation problem is `min t` s.t. `O X O ⪯ t² X`, `X = diag(x) ≻ 0`
   (§3). Its level sets in `x` are spectrahedra (linear matrix inequalities), so it is
   **quasi-convex** and solved by bisection with an SDP oracle — it is *not* non-convex in
   the harmful sense, and it is *not* a generalised eigenvalue problem in general. The
   generalised-eigenvalue structure appears only when the diagonal restriction on `X` is
   dropped (then the answer is `ρ(O)`, Stein/Lyapunov) or in the ∞-norm Perron case
   (Theorem 3.4), where the minimum is exactly `ρ(A)` for nonnegative `A`.
3. There is **no smaller enclosing eigenvalue problem** arising from signs (§5). The exact
   constant is the `n × n` symmetric eigenproblem; `ρ(|O|)` is a same-size Perron problem
   that is *cheaper to certify but never smaller*. Genuine size reductions come from rank
   or block structure, not from signs.
4. What signs *do* buy, sharply (§4): the chain
   `ρ(O) = ‖O‖₂ ≤ ρ(|O|) ≤ max_i Σ_j |O_ij|`, where the **middle** inequality is the
   Perron refinement of Gershgorin (magnitude-only, Collatz–Wielandt/Osborne) and the
   **first** is the genuine sign gain, with an exact equality criterion (Lemma 4.3).
   Because `ρ(O) = ‖O‖₂` is the sharp constant, all of these are upper bounds on `τ`
   that are *worse*, not better — the "gain" is that the magnitude route overestimates.

Nothing in §2–§5 is new mathematics; §7 says precisely what is classical and what (little)
is not.

All numeric values quoted below were computed independently (`drafts/check_sign_structure*.py`),
and the results marked **✓Lean** are machine-checked in `drafts/TauSignStructure.lean`
(axioms: `propext`, `Classical.choice`, `Quot.sound` only — no `sorry`, no custom axiom).

---

## 1. Notation

`O ∈ ℝ^{n×n}` symmetric; `S = diag(s)` with `s > 0`; `X = S² = diag(x)`, `x = s² > 0`;
`‖·‖₂` spectral norm; `ρ(·)` spectral radius; `|O|` the entrywise absolute value
(nonnegative, symmetric); `‖·‖_∞` the induced ∞-norm (max absolute row sum);
`λ_min`, `λ_max` the extreme eigenvalues of a symmetric matrix; `(a)⁺ = max(a,0)`,
`(a)⁻ = max(−a,0)`.

---

## 2. (1) The 2-norm cannot be reduced by a positive diagonal similarity

### Proposition 2.1 (no-go)
Let `O` be **normal** (in particular symmetric) and `S = diag(s)`, `s > 0`. Then

  `‖S O S⁻¹‖₂ ≥ ρ(S O S⁻¹) = ρ(O) = ‖O‖₂`,

with equality at `S = cI`. Hence `inf_{s>0} ‖S O S⁻¹‖₂ = ‖O‖₂`, and the infimum is attained.

*Proof.* `S O S⁻¹` is similar to `O`, so the two have the same characteristic polynomial and
hence the same spectrum; in particular `ρ(S O S⁻¹) = ρ(O)`. For any matrix `B` and any
eigenpair `Bv = λv` with `v ≠ 0`, `|λ| ‖v‖ = ‖Bv‖ ≤ ‖B‖₂ ‖v‖`, so `ρ(B) ≤ ‖B‖₂`. For
symmetric `O`, `‖O‖₂ = max_i |λ_i(O)| = ρ(O)`. Chaining,
`‖S O S⁻¹‖₂ ≥ ρ(S O S⁻¹) = ρ(O) = ‖O‖₂`. At `S = cI`, `S O S⁻¹ = O`. ∎

### Corollary 2.2 (`τ = ‖O‖₂` is optimal, hence sign-blind)
For symmetric `O`, `sup_{‖x‖=1} |Q_O(x)| = ‖O‖₂` (Rayleigh quotient; Horn–Johnson §4.2),
and the supremum is attained at an eigenvector for the eigenvalue of largest modulus.
Therefore any `τ`
satisfying the certificate hypothesis must satisfy `τ ≥ ‖O‖₂`, and `τ = ‖O‖₂` is achieved.
Sign information cannot lower this number: the sign pattern of `O` is already fully encoded
in its eigenvalues.

Numerical witness of the *sharpness* in dimension 2, **✓Lean**
(`two_by_two_scaling_cannot_help`): for `O = [[0,1],[1,0]]` and `S = diag(a,b)` one has
`S O S⁻¹ = [[0, a/b],[b/a, 0]]`, whose Rayleigh quotient at `x = (1,1)` is
`(a/b + b/a)/2 ≥ 1`, i.e. at least the unscaled value `1 = ‖O‖₂`, for **every** positive
`a,b`.

### 2.3 What the classical theorem actually says (and why "signed" is the wrong lever)
The textbook statement that diagonal scaling can make a norm collapse concerns
**non-normal** matrices, and it is a statement about *eigenvalue conditioning*
(Householder; Bauer–Fike 1960; Varah 1975 "On a generalization of Gershgorin's theorem";
Osborne 1960 "matrix balancing"):

* If `A = V Λ V⁻¹` is diagonalisable then `inf_{S invertible} ‖S A S⁻¹‖₂ = ρ(A)`,
  attained at `S = V⁻¹`; the quantity `inf_{S diagonal} ‖S A S⁻¹‖₂` (the *diagonal-scaling
  radius*, ≥ `ρ(A)`) is the relevant one for diagonal `S`.
* Standard illustration: `A = [[1,1000],[0,1]]` has `‖A‖₂ = 1000.0005` but
  `ρ(A) = 1` and `inf_{t>0} ‖diag(1,t) A diag(1,t)⁻¹‖₂ = 1` (computed: `1.0005` on a grid
  of `t`, limited only by the grid). Likewise `B = [[0,1000,0],[0,0,1000],[0.001,0,0]]`
  (nonnegative!) has `‖B‖₂ = 1000`, `ρ(B) = 10`, and `‖·‖_∞` drops from `1000` (row sums)
  to exactly `10` under the scaling `S = diag(1/v)` with `v` the Perron vector.

The mechanism is **non-normality**, not signs: for normal matrices the spectrum is
invariant under unitary similarity and `‖·‖₂ = ρ`, leaving no room. A symmetric `O` is
normal, so the collapse is impossible; a *signed symmetric* matrix gains nothing that an
unsigned one would not. Mixed signs are visible to `τ` only through `ρ(O)` itself.

---

## 3. (2) The exact optimisation problem, its convexity, and where eigenvalues appear

### 3.1 Reduction to a linear matrix inequality
For `S = diag(s) ≻ 0` and any `t ≥ 0`,

  `‖S O S⁻¹‖₂ ≤ t`
  `⟺ (S O S⁻¹)ᵀ(S O S⁻¹) ⪯ t² I`
  `⟺ S⁻¹ O S² O S⁻¹ ⪯ t² I`
  `⟺ O S² O ⪯ t² S²`  (congruence by `S ≻ 0`: `A ⪯ B ⟺ SAS ⪯ SBS` for `S ≻ 0`; Sylvester)

and with `X = S² = diag(x)`, `x > 0`:

  **(P)** `min t` subject to `O X O ⪯ t² X`, `X ∈ 𝒟_{++} := {diag(x) : x > 0}`.

Since `X` is diagonal, `O X O = Σ_k x_k (O e_k)(O e_k)ᵀ` is **linear** in `x`. Hence:

### Proposition 3.2 (quasi-convexity, not convexity, and not a GEP)
(a) For each fixed `t`, the feasible set `{x > 0 : O diag(x) O ⪯ t² diag(x)}` is the
intersection of the open positive orthant with a spectrahedron, hence **convex**; the sets
are nested decreasing in `t`. Therefore **(P) is quasi-convex** and is solved exactly by
bisection on `t` with an LMI (SDP) feasibility oracle at each step.
(b) The objective `s ↦ ‖diag(s) O diag(s)⁻¹‖₂` is **not convex** in `s` (numerically:
18000 random midpoint tests on `3 × 3` symmetric `O` produced 1605 violations of
`f((s+t)/2) ≤ (f(s)+f(t))/2`, the worst by a factor `34.2`). Nor is `(P)` a convex program
in `(t,x)`: the constraint `t² diag(x) − O diag(x) O ⪰ 0` is bilinear in `(t², x)`, so it is
not an LMI in those variables. So the true statement is: *not convex, but quasi-convex with
LMI level sets after fixing `t`* — the good structure sits in the variables `x = s²`, not in
`s` or in `(t,x)` jointly.
(c) `(P)` is **not** a generalised eigenvalue problem in general. A GEP-shaped
reformulation appears in exactly two places: (i) if the restriction `X ∈ 𝒟_{++}` is
dropped, see Proposition 3.3, where the answer is `ρ(A)`; (ii) when `n = 2`, where the LMI
is `2 × 2` and everything is explicit (this is the only "small" case, and it is small
because of the dimension, not because of the diagonal cone). In general `(P)` is the
*diagonal Lyapunov / D-stability* feasibility problem (Barker–Berman–Plemmons 1978,
"Positive diagonal solutions to the Lyapunov equations"), solved by SDP, not by an
eigensolver. Note in particular that a rank-deficient `X = diag(x)` with small support is
*not* a legitimate restriction: for `X` supported on `S ⊊ {1..n}` the inequality is
required to hold on all of `ℝⁿ` and hence forces the compression of `OXO` to `S^⊥` to be
negative semidefinite.

*Proof of (a).* `x ↦ O diag(x) O − t² diag(x)` is affine in `x`; `⪯ 0` is a closed convex
constraint (the PSD cone is convex and the map is affine); `{x > 0}` is convex. Nestedness:
`t₁ ≤ t₂` gives `t₁²diag(x) ⪯ t₂²diag(x)` on `x > 0`, so the `t₂`-set contains the
`t₁`-set. Quasi-convexity = convex level sets + nestedness. ∎

`(b)` is immediate from the bilinearity; the non-convexity of `s ↦ ‖diag(s)A diag(s)⁻¹‖`
is classical (it is why matrix balancing is done in logarithmic coordinates).

### Proposition 3.3 (the unconstrained-`X` version is exactly the spectral radius)
For any square `A` and `t > 0`:

  `∃ X ≻ 0 : Aᵀ X A ⪯ t² X  ⟺  ρ(A) ≤ t.`

*Proof.* (⇒) For an eigenpair `Av = λv`, `|λ|² v*Xv = v*AᵀXAv ≤ t² v*Xv`, and `v*Xv > 0`,
so `|λ| ≤ t`. (⇐) If `ρ(A) < t`, the Stein equation `AᵀXA − t²X = −I` has a positive
definite solution (Stein 1951; Lyapunov), which satisfies the constraint; for `ρ(A) = t`
take the limit or perturb `A`. For symmetric `O` one can take `X = I` after diagonalising,
since `O² ⪯ t²I ⟺ ρ(O) ≤ t`. ∎

So the "generalised eigenproblem" the question expects exists, but **only after deleting
the diagonal restriction** — and then it returns `ρ(A)`, which by Proposition 2.1 is what
the diagonal-restricted problem already attains when `O` is symmetric. For symmetric `O`
problem **(P)** therefore collapses: `t* = ‖O‖₂`.

### Theorem 3.4 (the ∞-norm case: here diagonal scaling genuinely beats Gershgorin)
Let `A ⪰ 0` (entrywise), `A` irreducible. Then

  `min_{s>0} ‖diag(s) A diag(s)⁻¹‖_∞ = ρ(A)`,

the minimum being attained at `s = 1/v` where `v ≻ 0` is the Perron vector of `A`. Since
`s = 1` gives the maximum absolute row sum `‖A‖_∞` (the Gershgorin bound), the optimisation
over the diagonal-similarity class **strictly** improves Gershgorin unless all row sums
coincide.

*Proof.* `‖diag(s)A diag(s)⁻¹‖_∞ = max_i s_i (A s⁻¹)_i`; writing `w = s⁻¹ > 0` this is
`max_i (Aw)_i/w_i`. **Lower bound:** for any `w > 0`, let `μ = min_i (Aw)_i/w_i`; then
`Aw ≥ μw`, and the Collatz–Wielandt argument (perturb `A` to `A + εJ > 0`, use the positive
left Perron vector of `A + εJ`, let `ε ↓ 0`) gives `ρ(A) ≥ μ`. **Upper bound:** at the
Perron vector, `Av = ρv` gives `(Av)_i/v_i = ρ` for all `i`, so the value is `ρ(A)`. ∎

Numerically (`B` above, irreducible and nonnegative): `ρ(B) = 10`, `‖B‖_∞ = 1000`,
`‖diag(1/v) B diag(1/v)⁻¹‖_∞ = 10.000000` — an exact match to three-plus digits.

**Moral for the question as posed.** The intuition "optimal diagonal scaling beats
Gershgorin" is *correct*, but the theorem lives in the **∞-norm for nonnegative matrices**
(Collatz 1942; Wielandt 1950; Varga, *Matrix Iterative Analysis*, Ch. 2; Osborne 1960). It
does not transplant to the 2-norm of a symmetric matrix, where Proposition 2.1 forbids it.

---

## 4. What sign structure does buy: the sharp hierarchy

Write `r_i = Σ_j |O_ij|` (absolute row sums) and recall `ρ(|O|) ≤ max_i r_i`
(Collatz–Wielandt with weights `1`; equivalently the Gershgorin/Schur bound). Then, for
symmetric `O`,

  `‖O‖₂ = ρ(O) ≤ ρ(|O|) ≤ max_i r_i`.  (★)

The first inequality is where signs matter.

### Lemma 4.1 (deficit identity; the pointwise form of (★))
**✓Lean** for the inequality `|Q_O(x)| ≤ Σ_{i,j}|O_ij||x_i||x_j|`; the exact identity below is
verified numerically to machine precision (`check_sign_structure3.py` family), not formalised.
For symmetric `O` and any `x`,

  `Q_O(x) = Σ_{i,j} |O_ij| |x_i| |x_j| − 2 Σ_i (O_ii)⁻ x_i²
            − 4 Σ_{i<j, O_ij>0} O_ij (x_i x_j)⁻ − 4 Σ_{i<j, O_ij<0} |O_ij| (x_i x_j)⁺`.

Consequently `|Q_O(x)| ≤ Σ_{i,j}|O_ij||x_i||x_j| ≤ max_i r_i · ‖x‖²`, and the deficit
vanishes identically iff `x_i x_j ≥ 0` on positive entries and `x_i x_j ≤ 0` on negative
entries (i.e. `sign(x)` reproduces the sign pattern of `O`, up to vertex switching).

*Proof.* Diagonal terms: `O_ii x_i² = |O_ii| x_i² − 2(O_ii)⁻ x_i²`. Off-diagonal pair
`{i,j}`: the form contributes `2O_ij x_i x_j`, the absolute-value form `2|O_ij||x_i||x_j|`,
so the deficit is `2|O_ij||x_i||x_j| − 2O_ij x_i x_j`. If `O_ij > 0` this equals
`2O_ij(|x_ix_j| − x_ix_j) = 4O_ij (x_ix_j)⁻`; if `O_ij < 0`,
`2|O_ij|(|x_ix_j| + x_ix_j) = 4|O_ij|(x_ix_j)⁺`. Summing gives the identity, each deficit
being nonnegative, whence the inequality. ∎

The Lean file proves the two consequences used downstream: `|Q_O(x)| ≤ Σ|O_ij||x_i||x_j|`
(`abs_quadform_le_absO`), the weighted/Gershgorin bound (`quadform_le_diag_scaled`,
`quadform_le_rowsum`), and the exact statement behind "negative entries are free on the
positive orthant" (`posOrthant_gain`, Lemma 5.1).

### Lemma 4.2 (`ρ(O) ≤ ρ(|O|)`)
Let `O` be any square matrix. Then `ρ(O) ≤ ρ(|O|)`.

*Proof.* Let `Ox = λx`, `|λ| = ρ(O)`, `x ≠ 0`. Componentwise,
`ρ(O)|x_i| = |λ x_i| = |Σ_j O_ij x_j| ≤ Σ_j |O_ij||x_j|`, i.e. `|O||x| ≥ ρ(O)|x|` with
`|x| ≥ 0`, `|x| ≠ 0`. For nonnegative `B` and `u ≥ 0`, `u ≠ 0`: `Bu ≥ μu` implies
`ρ(B) ≥ μ` (Collatz–Wielandt: with `J` the all-ones matrix and `ε > 0`, `(B+εJ)u > μu`;
take the positive left Perron vector `y` of `B + εJ`: `ρ(B+εJ) yᵀu = yᵀ(B+εJ)u > μ yᵀu`,
so `ρ(B+εJ) > μ`; then `ε ↓ 0`). Apply with `B = |O|`, `u = |x|`. ∎

### Lemma 4.3 (exact equality criterion)
Assume `|O|` is irreducible and `O` symmetric. Then
`ρ(O) = ρ(|O|)` **iff** there is an eigenvector `x` of `O` with
`Ox = ±ρ(O) x`, `|x| ≻ 0`, and **no row cancellation**, i.e.
`|O| |x| = ρ(O)|x|` (equivalently `|(Ox)_i| = (|O||x|)_i` for every `i`).

*Proof.* (⇐) `|x| ≻ 0` and `|O||x| = ρ(O)|x|` make `|x|` an eigenvector of the nonnegative
irreducible `|O|` for the eigenvalue `ρ(O) > 0`; a nonnegative irreducible matrix has a
positive eigenvector only for its Perron root, so `ρ(|O|) = ρ(O)`.
(⇒) Let `Ox = λx`, `|λ| = ρ(O) = ρ(|O|)`. From Lemma 4.2's chain,
`|O||x| ≥ ρ(O)|x|` with `|x| ≥ 0`, `|x| ≠ 0`, and `ρ(|O|) = ρ(O)`, so
`|O||x| ≥ ρ(|O|)|x|`. If some component were strict, the Collatz–Wielandt bound (applied to
`|O| + εJ` and the positive left Perron vector) would give `ρ(|O|) > ρ(O)`, a contradiction;
hence `|O||x| = ρ(O)|x|`, which with irreducibility gives `|x| ≻ 0` (Perron–Frobenius), and
every triangle inequality in the chain is an equality, i.e. no row cancellation. ∎

Verified numerically on 120 random signed `4 × 4` sign patterns: **0 mismatches** between
`ρ(|O|) − ρ(O) < 10⁻⁷` and the existence of a non-cancelling eigenvector at `±ρ`
(`drafts/check_sign_structure3.py`, part B′).

### Corollary 4.4 (signature similarity suffices, imbalance does not force a gap)
If `D` is a signature matrix (`D_ii = ±1`) with `D|O|D = O`, then `O` and `|O|` are similar
and `ρ(O) = ρ(|O|)`: no sign gain. **The converse is false.** For
`O = [[0,1,−1],[1,0,1],[−1,1,0]]` (one "frustrated" triangle, unbalanced) the characteristic
polynomial is `−(λ−1)²(λ+2)`, so `ρ(O) = 2 = ρ(|O|) = ρ(K₃)` even though no signature
matrix satisfies `D|O|D = O`. The correct equality criterion is Lemma 4.3, not balance.

### Corollary 4.5 (Z-matrices: the magnitude bound is exact)
If `O = −B` with `B ⪰ 0` (all off-diagonal entries of `O` nonpositive — the *opposite* of
"mixed sign"), then `|O| = B`, `ρ(O) = ρ(B) = ρ(|O|)`, so `‖O‖₂ = ρ(|O|)`; and
`λ_min(D + O) = λ_min(D − B)` is controlled exactly by the M-matrix criterion
`ρ(D⁻¹B) < 1` (Perron–Frobenius; Varga). In that regime the *only* available refinement is
Gershgorin → Perron (Theorem 3.4), i.e. the magnitude-only step of (★).

### Where the gaps are, numerically
| `O` | `ρ(O) = ‖O‖₂` | `ρ(|O|)` | `max_i r_i` |
|---|---|---|---|
| signed `C₄` with one negative edge | `√2 ≈ 1.4142` | `2` | `2` |
| frustrated triangle `[[0,1,−1],[1,0,1],[−1,1,0]]` | `2` | `2` | `2` |
| `[[2,3,−1],[3,1,4],[−1,4,0.5]]` | `5.5890` | `6.6155` | `8` |
| signed `C₄`, but `O = −B`, `B ⪰ 0` (Z-matrix) | `ρ(B)` | `ρ(B)` | `max_i r_i` |

Both inequality steps in (★) are thus separately visible, and every entry of the table is
`≥ ‖O‖₂`: **the sign-aware refinements only ever move `τ` in the wrong direction** — they
shrink an *overestimate*, they do not go below the sharp constant.

---

## 5. (3) Is there a smaller enclosing eigenvalue problem?

**(a) From signs alone: no.** The exact constant is `sup_{‖x‖=1}|Q_O(x)| = ρ(O)`, an `n × n`
symmetric eigenproblem. A sign pattern is not a sparsity structure: it constrains only the
*signs* of the entries, not their sizes, so it neither decouples the matrix nor shrinks it.
Precisely: for a fixed sign pattern the matrices realising it form a cone `𝒞` on which
`sup_{O∈𝒞, |O_ij| fixed} ρ(O)` depends on the actual signed spectrum (Lemma 4.3), and no
functional of the magnitudes alone can equal it — the magnitude-only functionals are
exactly the upper bounds in (★), each of them attained only for sign patterns that are
switching-consistent with a Perron vector.

**(b) `ρ(|O|)` is same-size and weaker.** It is a Perron eigenproblem — computable by
monotone power iteration with no cancellation, hence a cheaper *certificate* — but by (★)
it is `≥ ρ(O)`. It cannot be the "smaller enclosing eigenvalue problem" asked for; it
encloses `ρ(O)` from above at the same dimension.

**(c) Genuine dimension reductions exist, and they are not sign-based.** If `O = ABᵀ` has
rank `r ≪ n` then `‖O‖₂` is the root of an `r × r` eigenproblem (`O² = ABᵀBA`, or the thin
SVD); if `O` is block/banded up to permutation the eigenproblem splits. Sign patterns often
*co-occur* with low rank (`O_ij = u_i v_j`-type structure), which is presumably why signs
feel informative, but the reduction is bought by rank, not by signs.

**(d) The one exact sign-based simplification: a positive ground state.** If `O` has
nonpositive off-diagonal entries, Perron–Frobenius gives a positive ground vector of
`D + O` (equivalently `D − B` with `B ⪰ 0` is a Z-matrix), and the eigenproblem becomes a
Perron problem on the cone. The relevant formal fact is

### Lemma 5.1 (negative entries are free on the positive orthant) **✓Lean**
If `x_i x_j ≥ 0` for all `i,j` (in particular if `x ≥ 0`) then
`Q_O(x) ≤ Σ_{i,j} (O_ij)⁺ x_i x_j`: the negative part of `O` can be deleted.

*Proof.* Termwise, `O_ij x_i x_j ≤ max(O_ij,0) x_i x_j`: if `O_ij ≥ 0` it is an equality; if
`O_ij < 0` the left side is `≤ 0 ≤` the right side. ∎

Using 5.1 inside the certificate requires the *extra hypothesis* that the minimising vector
is sign-definite; that hypothesis, not the sign pattern per se, is what discards the
negative entries — and it is exactly why mixed signs are unhelpful: a mixed-sign `O` has
ground vectors of both signs in general, while a Z-matrix has a positive one.

**(e) A different (sign-free) scaling remark about the certificate.** Congruence by `D^{-1/2}`
turns positivity into an exact criterion: with `K = D^{-1/2} O D^{-1/2}`,
`D + O ≻ 0 ⟺ I + K ≻ 0 ⟺ λ_min(K) > −1` (Sylvester's law of inertia; checked numerically on
200 random cases, 0 mismatches). This criterion is strictly stronger than the sufficient
condition `d > ‖O‖₂` (which gives `λ_min(K) ≥ −‖O‖₂/d > −1`), but it is a *criterion*, not a
numeric bound: the congruence alone yields the numeric lower bound
`λ_min(D+O) ≥ d(1 − ‖K‖₂)` **only under the extra hypothesis `‖K‖₂ ≤ 1`** (proof: `I + K ⪰
(1−‖K‖₂)I`, so `D + O = D^{1/2}(I+K)D^{1/2} ⪰ (1−‖K‖₂)D ⪰ d(1−‖K‖₂)I`, where the last step
needs the scalar `1−‖K‖₂ ≥ 0`; verified on 231 random cases with `‖K‖₂ < 1`). Without that
hypothesis the analogous inequality is **false** (`d(1−‖K‖₂)` can exceed `λ_min(D+O)` when
`1−‖K‖₂ < 0` — checked and rejected numerically). No sign information is used anywhere here.
When the hypothesis fails one falls back on Weyl: `λ_min(D+O) ≥ d − ‖O‖₂`.

---

## 6. What is machine-checked, and what is not

`drafts/TauSignStructure.lean` (compiles with Lean 4.19.0 + mathlib v4.19.0; axiom audit
shows only `propext`, `Classical.choice`, `Quot.sound`):

| theorem | content |
|---|---|
| `abs_quadform_le_absO` | `|Q_O(x)| ≤ Σ|O_ij||x_i||x_j|` (Lemma 4.1's inequality) |
| `two_mul_abs_le_weighted` | the weighted AM–GM step `2ab ≤ (v_j/v_i)a² + (v_i/v_j)b²` |
| `quadform_le_diag_scaled` | the diagonal-similarity / weighted Collatz–Wielandt bound: `Σ_j|O_ij|v_j ≤ τ v_i ∀i ⟹ |Q_O(x)| ≤ τ‖x‖²` |
| `quadform_le_rowsum` | its `v ≡ 1` case: the Gershgorin/Schur row-sum bound |
| `posOrthant_gain` | Lemma 5.1 |
| `two_by_two_scaling_cannot_help` | Corollary 2.2's witness: `1 ≤ (a/b + b/a)/2` |

The `quadform_le_diag_scaled` hypothesis is exactly `‖diag(v)|O|diag(v)⁻¹‖_∞ ≤ τ`; the
optimal weight is the Perron vector of `|O|` (giving `ρ(|O|)`), the trivial weight `v ≡ 1`
gives the row-sum bound. The Perron optimality itself (Theorem 3.4) is *not* formalised.

**Not formalised:** Proposition 2.1 and Corollary 2.2 in full (they need the spectral
theorem / min–max characterisation of `‖·‖₂`; the available mathlib build cache in this
workspace lacks the spectral library). It is a five-line textbook argument (given above) and
its conclusion is confirmed numerically (20 000-point random search plus a coordinate-descent
minimiser over `log s`: the minimum found is `‖O‖₂` to machine precision, gain `0.0`).

---

## 7. Novelty assessment

Nothing above is new mathematics. Attributions:

* Proposition 2.1 / the diagonal-scaling radius: classical (Householder; Bauer–Fike 1960;
  Varah 1975; Osborne 1960). The ∞-norm/Perron theorem 3.4: Collatz (1942), Wielandt (1950),
  Varga's *Matrix Iterative Analysis*.
* Reduction of the scaling problem to an LMI and its quasi-convexity: standard convex
  optimisation (Boyd–Vandenberghe, *Convex Optimization*, §4 and the geometric-programming
  chapter); the diagonal-Lyapunov problem is Barker–Berman–Plemmons (1978).
* `ρ(O) ≤ ρ(|O|)`: classical (Berman–Plemmons, *Nonnegative Matrices in the Mathematical
  Sciences*; Horn–Johnson §8). The equality criterion Lemma 4.3 is folklore in signed-graph
  spectral theory (Harary 1953 for balance; Zaslavsky for signed graphs); the precise
  "no row cancellation at a Perron vector" form is what the chain argument yields and is
  stated here for completeness.
* Weyl (1912) for the two-sided eigenvalue bound; Ostrowski/Brauer for the sign-refined disc
  theorems (still magnitude-only, hence still `≥ ρ(O)` and not useful here).

What is specific to this project is only a **negative** result and a bookkeeping service:
within the `λ_min(D+O) ≥ d − τ` formulation, `τ` cannot be improved by exploiting the signs
of `O_ij`, because the sharp constant `‖O‖₂` is attained; every sign-aware surrogate in the
literature is an upper bound on `‖O‖₂` and hence a *weaker* `τ`. Consequently the only
routes to a better certificate are (i) a better lower bound on the diagonal part, (ii) a
larger `d` obtained by absorbing part of `O` into a rescaled diagonal, or (iii) replacing the
uniform statement by a cone-restricted one (Lemma 5.1), which requires a Perron–Frobenius
hypothesis on the signs rather than merely "mixed signs".

---

### Reproduce

```
cd /mnt/e/CalcChemistry
export PATH=$PWD/.tools/lean-4.19.0/bin:$PATH
lake env lean drafts/TauSignStructure.lean     # exits 0, no warnings
python3 drafts/check_sign_structure.py         # numeric checks for §2, §4
python3 drafts/check_sign_structure2.py        # balance vs gap, Perron scaling
python3 drafts/check_sign_structure3.py        # equality criterion, ∞-norm theorem
python3 drafts/check_convexity.py             # §3: non-convexity in s, no 2-norm gain
```