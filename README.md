# CalcChemistry: verified theory and research audit

The previous status report overstated unverified results. This repository is
being rebuilt around finite-dimensional energy bounds and their conditional
propagation to chemical reaction rates. See STATUS_REPORT.md and docs/.

No universal polynomial-time molecular solver or validated ALD predictor has
been established. Historical surveys and scripts are exploratory material,
not evidence of a theorem or a chemical measurement.

The user-supplied targets are Al2O3/MoCl5 and CaF2 crystal deposition. Their exact meaning, surface
model, precursor/co-reactant pair, and operating conditions remain to be
specified. The user confirmed both spellings directly (MoCl5 on request, and
CaF2 crystal deposition); earlier drafts of this repository carried `CaF` and
`MoCl4`, which were transcription errors on the user's side, now corrected.

The working manuscript is `paper_main.tex` with five included sections and
`references.bib` (eight DOI records checked against Crossref). It gives full
written proofs under explicit finite-model assumptions; it is not a claim of
a new universal solver or a completed ALD calculation.

## Reproduction

Python 3.10 was used for the numerical checks. Create a local environment and
install `requirements-verified.txt` for NumPy, SciPy and Matplotlib. Then run:

    python -B -m unittest discover -s tests -v
    python -B richardson_h2_verification.py
    python -B -m experiments.certificates

Lean and mathlib are pinned to v4.19.0. For ordinary network access:

    lake update
    lake exe cache get
    lake build
    lake env lean scripts/CheckAxioms.lean

If Git smart-HTTP fails after the mathlib checkout succeeds, use
`scripts/bootstrap_archives.py` to obtain the eight upstream-pinned source
archives. It verifies existing source files before reuse, records archive
hashes, and writes local dependency paths into the root manifest. This is a
transport fallback, not a change to mathlib's source or proof checking.
The root `lake-manifest.json` in this workspace records that fallback;
a fresh checkout needs the archive bootstrap before building with it.

The manuscript can be compiled with Tectonic 0.15.0 or a normal LaTeX/BibTeX
installation. The complete local gate is:

    python scripts/verify.py

Build artifacts are written to ignored `build/`. PDF generation is only a
compile check; the requested deliverable remains the LaTeX source. Actual run
status and theorem coverage are in `STATUS_REPORT.md` and
`docs/VERIFICATION_SCOPE.md`. Julia is corrected source only and remains
untested until a Julia runtime is used.
