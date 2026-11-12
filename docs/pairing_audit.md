# Pairing and finite-certificate numerical audit

## Scope and corrected pairing model

The historical `richardson_h2_*` filenames describe a two-orbital reduced pairing toy model. They do not contain molecular integrals or an ab initio H₂ error certificate. `richardson_h2_proof.md` now derives the complete six-state spectrum directly.

| Issue | Correction |
| --- | --- |
| Pair creation checked before annihilation | Apply annihilation first so the `p=q` self-pairing term contributes `-g` on each paired determinant. |
| Missing diagonal energy | The paired `H_z` block is `[[eps1-eps2-g, -g], [-g, eps2-eps1-g]]`. |
| Mixed energy origins | `H_N = H_z + (eps1+eps2) I`; Richardson roots are number-convention energies. |
| Negative secular sign | For interaction `-g P†P`, `1 = g sum_p 1/(2 eps_p - E_pair)` away from poles. |
| Incomplete state count | Two paired states plus four unpaired states give all six determinants. |
| Singlet label too broad | The paired block is the seniority-zero subspace, not the whole singlet space. |
| Degeneracies/poles | Use matrix residuals and eigenspace projection; a missing rational residual is not evidence of success. |
| Molecular and integrability overclaims | Remove claims that pairing dominance makes the reduction exact, and remove unsupported commuting-charge/completeness formulas. |

For `eps1=-1.0`, `eps2=-0.5`, `g=0.3`, the previous `H_z` value **-0.583095189484530** is replaced by **-0.883095189484530**. Its corresponding `H_N` value is **-2.383095189484530**. Omitting self-pairing changes the operator, not merely its common energy origin, because unpaired states are unaffected.

The actual Python CLI run passed ten parameter cases. For this reference case it returned energy and full-spectrum errors `2.220e-16`, maximum paired-vector residual `1.570e-16`, and normalization error `0.0`. These are observed floating-point errors, not a universal roundoff guarantee. The pairing Python implementation and its tests were not edited during this certificate/documentation pass.

## Block-energy enclosure and failure modes

For a finite Hermitian matrix

\[
H=\begin{pmatrix}A&B\\B^*&C\end{pmatrix},\quad
 a=\lambda_{\min}(A),\quad c=\lambda_{\min}(C),\quad b=\|B\|_2,
\]

write a normalized vector as `(x,y)`. The quadratic form satisfies

\[
\langle (x,y),H(x,y)\rangle
\ge a\|x\|^2+c\|y\|^2-2b\|x\|\|y\|.
\]

Minimizing the right side is a two-dimensional eigenvalue problem, while a lowest active-block vector is a variational trial state. Therefore

\[
L=\frac{a+c-\sqrt{(c-a)^2+4b^2}}2
\le\lambda_{\min}(H)\le U=a.
\]

The lower bound is sharp for scalar real/complex blocks. If `c>a`, rationalizing the square root yields

\[
0\le a-\lambda_{\min}(H)\le a-L
=\frac{2b^2}{\sqrt{(c-a)^2+4b^2}+c-a}
\le\frac{b^2}{c-a}.
\]

The positive-gap estimate does not apply when `c<=a`. For `H=diag(1,-2)` with the first state active, the active residual is exactly zero, but the ground energy is `-2`. The implementation returns `[-2,1]` and `gap_error_bound=None`, preventing a false ground-state conclusion from an excited eigenvector's zero residual.

The implementation uses a scaled rationalized correction to avoid squaring very small/large couplings first. Tests reproduced the old failures at `b=1e-200` and `b=1e200` before the fix. Invalid matrix shapes, nonfinite entries, non-Hermitian inputs, empty splits, and Boolean dimensions are rejected. Hermiticity is checked at absolute tolerance `1e-12`; it is a numerical input check, not an exact arithmetic proof.

## Barrier and conditional rate intervals

Given energy enclosures `E_R in [L_R,U_R]` and `E_T in [L_T,U_T]`, subtraction gives

\[
\Delta E\in[L_T-U_R,U_T-L_R].
\]

For a fixed nonnegative prefactor `A`, fixed positive thermal energy `theta`, and the assumed rate law `k=A exp(-Delta E/theta)`, monotonicity gives

\[
k\in[Ae^{-U/\theta},Ae^{-L/\theta}].
\]

`barrier_interval(reactant_interval, transition_interval)` preserves negative barriers. `rate_interval(barrier, prefactor, thermal_energy)` reverses the endpoints. Both require ordered finite real intervals; the rate helper rejects negative prefactors, nonpositive thermal energy, nonfinite inputs, and positive rates that overflow or underflow to zero. Zero prefactor is allowed. The input barrier and `theta` must share units: `k_B T` for per-particle energies or `R T` for molar energies. Prefactor and temperature uncertainty, entropy corrections, and the validity of a reaction model are not certified here.

The synthetic example has reactant interval `[-100.1,-99.9] kJ/mol`, transition interval `[-80.2,-79.8] kJ/mol`, `T=300 K`, and `A=1e6 s^-1`. It produces a barrier approximately `[19.7,20.3] kJ/mol` and rates approximately `[292.09561581036974,371.52862369194713] s^-1`. These inputs are invented test parameters, explicitly marked `synthetic_only`; they are not molecular measurements or computed chemistry data. Endpoint comparisons show a maximum floating-point rate violation `1.7053025658242404e-13`, below the recorded tolerance. It is deliberately not rounded away or reported as a rigorous enclosure.

## Fixed finite-state CTMC check

For row generators `Q` and `Q_tilde`, with nonnegative off-diagonal entries and zero row sums, set

\[
\eta=\tfrac12\max_i\sum_j|Q_{ij}-\widetilde Q_{ij}|.
\]

For the same initial probability distribution `p0`, Duhamel's identity gives

\[
p_0(e^{tQ}-e^{t\widetilde Q})
=\int_0^t p_0 e^{sQ}(Q-\widetilde Q)e^{(t-s)\widetilde Q}\,ds.
\]

Right multiplication by a stochastic matrix contracts the `l1` norm, and a convex combination of generator-difference rows has `l1` norm at most `2 eta`. Thus

\[
\operatorname{TV}(p_0e^{tQ},p_0e^{t\widetilde Q})\le t\eta.
\]

The experiment uses two explicitly stored three-state generators, rates in `s^-1`, times in seconds, the three point-mass initial distributions and `[0.2,0.3,0.5]`, at times `[0,0.01,0.1,1,5]`. `scipy.linalg.expm` computes all 20 comparisons. Tests check generator validity, probability conservation, nonnegativity, independently recomputed TV, and the bound. This uses fixed generators only; it is not a numerical test of a time-dependent or molecular kinetic model.

## Executed evidence and reproduction

```sh
/usr/bin/python3 -B -m unittest discover -s tests -v
/usr/bin/python3 -B richardson_h2_verification.py
/usr/bin/python3 -B -m experiments.certificates
```

Observed results:

- **30 tests passed in the final shared-workspace run:** 15 certificate/report tests, 12 unchanged pairing tests, and 3 audit-counterexample tests added independently. The final run enabled Python warnings as errors and NumPy `seterr(all='raise')`; all passed. The earlier certificate/pairing-only run passed 27 tests.
- **360 random block splits:** real symmetric and complex Hermitian matrices of dimensions 2 through 10, four samples of each kind per dimension, every nontrivial leading split, seed `20261112`, NumPy `PCG64`.
- Maximum random enclosure violation: **`4.440892098500626e-16`**; maximum applicable positive-gap violation: **`0.0`**.
- Nine sharp two-by-two cases; maximum sharpness error **`8.881784197001252e-16`**.
- **20 CTMC comparisons:** maximum `TV - t eta`, clipped below by zero, **`0.0`**.
- JSON generated by actual execution with Python `3.10.12`, NumPy `2.2.6`, SciPy `1.15.3` and overall `passed=true`.

`experiments/certificate_results.json` stores every random split's dimension, sample, scalar block data, numerical eigenvalue reference, tolerance and errors; it also stores the fixed example parameters and per-case CTMC probabilities. Random matrices can be reconstructed from the recorded generator, seed, construction rule, sample order (real then complex, increasing dimension), and library version. `exact` is the legacy field name for a floating-point dense diagonalization reference, not an exact number. Acceptance uses `5e-12 * max(1, ||H||_2)` for block tests and `5e-13` for CTMC checks. The script enables NumPy floating-point exceptions and refuses nonfinite JSON output.

## Julia correction status

**UNTESTED: no Julia executable was found.** Source changes include the `p=q` term, positive secular sign, number-to-quasi-spin energy conversion, stable paired vectors with explicit `g`, full spectrum and both vector residual checks, and degenerate eigenspace comparison. The script's success message is conditional on assertions. The optional BigFloat check uses 256-bit scalar polynomial/vector algebra; it does not claim a generic BigFloat matrix eigensolver, symbolic verification, or rigorous interval arithmetic.

Once Julia and `Combinatorics` are available, run:

```sh
julia richardson_h2_verification.jl
julia richardson_h2_verification.jl --bigfloat
```

Neither command was executed here; no Julia result is included in the passed-test count. All claimed validation remains finite, synthetic, and floating point.
