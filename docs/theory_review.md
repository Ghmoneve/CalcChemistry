# Theory review: what is proved, what is conditional, and what remains open

## 1. Scope

The project cannot currently claim a universal polynomial-time solver for molecular electronic structure. Fixed-basis electronic structure and N-representability have QMA-complete worst-case promise-problem formulations (sources OGorman2022 and Liu2007 in `references.bib`; evidence in `docs/reference_audit.md`). This is a complexity-theoretic barrier, not an unconditional proof that polynomial algorithms cannot exist: no unproved separation such as BQP ≠ QMA is asserted here.

The defensible contribution is a certificate framework: if a finite Hamiltonian is split into an active block and a complementary block, computable spectral information can produce an energy interval. Interval subtraction gives a reaction-energy or barrier interval, and a stated transition-state formula can propagate that interval to a conditional rate interval.

## 2. Block energy enclosure

Let

\[
H=\begin{pmatrix}A&B\\B^*&C\end{pmatrix}
\]

be Hermitian, with
\(a=\lambda_{\min}(A)\), \(c=\lambda_{\min}(C)\), and \(b=\|B\|_2\). For every unit vector \((x,y)\),

\[
\langle(x,y),H(x,y)\rangle
\ge a\|x\|^2+c\|y\|^2-2b\|x\|\|y\|.
\]

The right-hand side is the quadratic form of

\[
M=\begin{pmatrix}a&-b\\-b&c\end{pmatrix}
\]

on \((\|x\|,\|y\|)\). Therefore

\[
\ell:=\frac{a+c-\sqrt{(a-c)^2+4b^2}}2
\le \lambda_{\min}(H).
\]

The variational principle, using vectors supported in the active block, gives
\(\lambda_{\min}(H)\le a\). Thus

\[
\boxed{\ell\le E_0(H)\le a.}
\]

If \(c>a\), then \(\ell\ge a-b^2/(c-a)\), so

\[
0\le a-E_0(H)\le \frac{b^2}{c-a}.
\]

This is a finite-dimensional theorem. It is not a claim that \(c\) or \(b\) can be obtained cheaply for an arbitrary molecular Hamiltonian.

## 3. Reaction energies and barriers

If \(E_{TS}\in[L_{TS},U_{TS}]\) and \(E_R\in[L_R,U_R]\), then

\[
L_{TS}-U_R\le E_{TS}-E_R\le U_{TS}-L_R.
\]

The same statement applies to a free-energy model after all zero-point, thermal, solvation, and entropic terms have themselves been bounded. Electronic energy alone is not a free-energy or rate certificate.

## 4. Conditional rate propagation

Under a separately stated Eyring/TST model with fixed prefactor \(K\ge0\), gas constant \(R>0\), temperature \(T>0\), and barrier \(B\in[L,U]\),

\[
K e^{-U/(RT)}\le K e^{-B/(RT)}\le K e^{-L/(RT)}.
\]

If the prefactor is uncertain, it needs an interval too; multiplying an electronic-energy interval by a guessed prefactor is not a rigorous rate prediction. The current Lean source formalizes interval subtraction and real exponential propagation with an interval prefactor in `CalcChemistry/CertifiedTheory.lean`. Compilation and axiom checks, rather than the presence of theorem declarations, determine verification status; see `STATUS_REPORT.md`.

## 5. Application boundary for ALD

For Al2O3/MoCl5, a future calculation must specify: surface polymorph and facet, hydroxyl termination, slab thickness and boundary conditions, MoCl5 charge and spin state, co-reactant, coverage, temperature, and a reproducible electronic-structure Hamiltonian. The user has confirmed the formula is MoCl5 and that the target is CaF2 crystal deposition; the earlier `MoCl4` and `CaF` were transcription errors on the user's side.

The certificate can then compare competing pathways only if each pathway receives a sound interval. It cannot certify a pathway from a hand-written geometry, a hard-coded DFT reference, or a guessed HOMO-LUMO gap.

## 6. Explicit limitations

1. The block theorem is an enclosure, not a polynomial-time algorithm.
2. A finite-basis certificate is not a basis-set-completeness certificate.
3. A variational upper bound alone does not give an error bar; a lower bound or residual theorem is required.
4. The previous repository SOS tables were hard-coded and are removed from the evidence base.
5. The previous persistent-homology example used a hard-coded reference count and is not independent DFT validation.
6. Born--Oppenheimer asymptotic results require gap, regularity, scaling, and constants; a numerical `epsilon^5` calculation alone is not an error bound.

## 7. Research target

The concrete research problem is now: construct and validate an active-space certificate for one specified ALD elementary step, then test whether the interval is narrow enough to rank two pathways. A negative result is scientifically useful: it identifies which missing information, rather than a generic exponential-wall slogan, controls the calculation.
