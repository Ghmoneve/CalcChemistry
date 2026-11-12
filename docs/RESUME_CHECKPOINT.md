# Checkpoint and resume record

The stopped jobs were restored by purpose, not by replaying obsolete commands.

- Lean 4.19.0 download: resumed with HTTP range, completed, extracted, version checked; linked into elan under the pinned toolchain name.
- Tectonic 0.15.0 download: resumed with HTTP range, completed, extracted, version checked.
- Initial mathlib setup: relaunched after the missing toolchain was fixed. Uses mathlib v4.19.0, source revision c44e0c8ee63ca166450922a373c7409c5d26b00b. Git smart-HTTP failed for plausible; `scripts/bootstrap_archives.py` restored all eight transitive dependencies from official commit-specific codeload archives, verified archive bytes on restart, and recorded URL/revision/SHA256 in `build/dependency-archives.json`. Root manifest uses local path entries; upstream mathlib sources and its lockfile remain unchanged.
- Official cache executable built successfully. Targeted cache arguments use module names when invoked from the project root. The optional fallback `scripts/CacheSelected.lean` calls upstream cache hashing/download/unpack functions without requiring a Git tag for optional ProofWidgets JS; it does not alter proof sources or kernel checking.
- Interrupted source verification: relaunched as a bounded eight-source metadata/evidence task.
- Interrupted full-theory manuscript: relaunched as a proof-complete finite-model manuscript task.
- Interrupted numerical branch: resumed from the corrected pairing implementation, with independent CAR tests and added certificate tests.
- Completed numeric audit: preserved, not rerun as brainstorming.

No prior research result is inferred from a background job's successful exit alone.
See STATUS_REPORT.md and build logs for final verification status.
