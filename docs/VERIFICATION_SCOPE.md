# Verification boundary

## Authoritative entry points for this revision

- `paper_main.tex` and its included sections: full written arguments within stated finite-model hypotheses; manuscript revision in progress.
- `CalcChemistry.lean`: imports only `CalcChemistry/Basic.lean`, `CalcChemistry/BornOppenheimer.lean`, and `CalcChemistry/CertifiedTheory.lean`.
- `scripts/CheckAxioms.lean`: prints theorem dependencies; `sorryAx` or additional physical axioms are not acceptable in the verified closure.
- `scripts/verify.py`: runs the numerical suite, pairing CLI, Lake build, axiom report and LaTeX build with real exit-code checking.
- `docs/audit_numerics.md`: evidence behind withdrawals of old numerical claims.

## Important distinctions

1. A successful build is not proof that the formal model describes a molecule.
2. A theorem over integers or natural numbers does not formalize a real-energy statement. The final core uses mathlib real numbers and exponential functions.
3. A finite-basis theorem does not bound basis-set, geometry, relativistic, thermal or surface-model errors unless each is independently enclosed.
4. A float regression test is not rigorous outward-rounded interval arithmetic.
5. A solved reduced pairing model is not an ab initio H2 or ALD calculation.
6. A conditional rate bound does not prove transition-state theory, discover a transition state, or enumerate all reactions.
7. Historical root-level Richardson/Yangian Lean files are archived, excluded drafts, not part of the proof closure.
8. Old SOS, topology and BO scripts contain assumed/constructed data. Their prior output is not an independent chemistry benchmark.

## Intended formal coverage

The core source proves real interval subtraction and additive-error propagation;
a scalar positive-semidefinite-form lemma, the sharp two-block scalar bound,
a completion-of-squares lower bound and its normalized version;
transfer from a uniformly bounded Rayleigh objective with an attained minimum;
monotone propagation through a positive-temperature real exponential with an
interval prefactor; relaxation order; disjoint-interval ranking; and the correct
paired-block characteristic polynomial/root identity.

The matrix-to-quadratic-form operator inequality, the full finite-dimensional
spectral theorem application, the continuous nuclear problem, BO asymptotics,
TST physics, and CTMC semigroup bound are not claimed as fully formalized by
these scalar lemmas. Check `STATUS_REPORT.md` for the actual completed runs.
