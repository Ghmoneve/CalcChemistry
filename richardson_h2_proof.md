# Two-level reduced pairing model: exact finite-sector derivation

The historical filename contains `h2`; the model here is a **synthetic reduced pairing Hamiltonian**, not an ab initio H₂ calculation. No molecular integrals, nuclear repulsion, geometry, basis-set limit, or molecular error bound are supplied. Exact algebra for this model does not establish chemical accuracy. Numerical comparisons below use ordinary floating-point arithmetic, not rigorous interval arithmetic.

## 1. Operators and energy conventions

Order the four fermionic modes as \((1\uparrow,1\downarrow,2\uparrow,2\downarrow)\), with canonical anticommutation relations. Define

\[
P_p^\dagger=a_{p\uparrow}^\dagger a_{p\downarrow}^\dagger,
\quad P_p=a_{p\downarrow}a_{p\uparrow},
\quad S_p^z=\tfrac12(n_{p\uparrow}+n_{p\downarrow}-1).
\]

The CAR give \([S_p^z,P_q^\dagger]=\delta_{pq}P_p^\dagger\), \([S_p^z,P_q]=-\delta_{pq}P_p\), and \([P_p^\dagger,P_q]=2\delta_{pq}S_p^z\). In particular,

\[
P_p^\dagger P_p=n_{p\uparrow}n_{p\downarrow}.
\]

For finite real orbital parameters and attractive coupling \(g\ge0\), distinguish

\[
H_z=\sum_{p=1}^2 2\varepsilon_p S_p^z
       -g\sum_{p,q=1}^2P_p^\dagger P_q,
\qquad
H_N=\sum_{p=1}^2\varepsilon_p(n_{p\uparrow}+n_{p\downarrow})
       -g\sum_{p,q=1}^2P_p^\dagger P_q.
\]

Both sums include **\(p=q\)**. The exact operator identity is

\[
H_N=H_z+(\varepsilon_1+\varepsilon_2)I.
\]

Thus a Richardson pair energy \(E_{\rm pair}\) in the number convention corresponds to \(E_z=E_{\rm pair}-(\varepsilon_1+\varepsilon_2)\). In particular, the empty vacuum has energy zero under \(H_N\), but energy \(-\varepsilon_1-\varepsilon_2\) under \(H_z\). A common orbital offset changes \(H_N\) by twice the offset in the two-electron sector while leaving \(H_z\) unchanged there.

## 2. Paired block, self-pairing, and complete spectrum

The two-electron space has dimension \(\binom42=6\). Its paired determinants are

\[
|1\bar1\rangle=P_1^\dagger|0\rangle=|1100\rangle,
\qquad |2\bar2\rangle=P_2^\dagger|0\rangle=|0011\rangle.
\]

Annihilating a pair and then creating it at the same orbital gives the same determinant. Consequently each occupied pair receives \(-g\) from \(p=q\). Checking whether the target orbital is empty **before** annihilation incorrectly removes this contribution. With \(d=\varepsilon_1-\varepsilon_2\),

\[
H_z^{\rm pair}=\begin{pmatrix}d-g&-g\\-g&-d-g\end{pmatrix},
\qquad
H_N^{\rm pair}=\begin{pmatrix}2\varepsilon_1-g&-g\\-g&2\varepsilon_2-g\end{pmatrix}.
\]

Therefore

\[
\det(H_z^{\rm pair}-\lambda I)=(\lambda+g)^2-d^2-g^2,
\quad
\boxed{\lambda_\pm=-g\pm\sqrt{d^2+g^2}},
\]

\[
\boxed{E_\pm=\varepsilon_1+\varepsilon_2-g\pm\sqrt{d^2+g^2}}.
\]

The **four** unpaired determinants are

\[
|1010\rangle,\ |1001\rangle,\ |0110\rangle,\ |0101\rangle.
\]

Each has one electron at each spatial orbital. All pair annihilators kill it and both \(S_p^z\) vanish. Their energies are zero under \(H_z\), and \(\varepsilon_1+\varepsilon_2\) under \(H_N\). This accounts for every state:

\[
\operatorname{spec}(H_z|_{N=2})=\{\lambda_-,0,0,0,0,\lambda_+\}.
\]

For \(g\ge0\), \(\lambda_-\le0\le\lambda_+\). When \(d\ne0\), the upper paired root is above the four zero-energy states: it is **not** the first excited energy of the full sector. The paired block is the seniority-zero subspace, not the whole spin-singlet subspace; the latter also contains a singly occupied singlet.

## 3. One-pair secular equation and eigenvectors

Write a paired state as \(|\psi\rangle=c_1|1\bar1\rangle+c_2|2\bar2\rangle\). The number-convention eigenvalue equation reads

\[
(2\varepsilon_p-E)c_p=g(c_1+c_2),\qquad p=1,2.
\]

If \(g>0\), \(c_1+c_2\ne0\), and neither denominator vanishes, then

\[
c_p\propto\frac1{2\varepsilon_p-E},
\qquad
\boxed{1=g\left(\frac1{2\varepsilon_1-E}+\frac1{2\varepsilon_2-E}\right)}.
\]

The sign is **positive** for the attractive interaction \(-gP^\dagger P\). Clearing denominators away from the poles gives

\[
(2\varepsilon_1-E)(2\varepsilon_2-E)
-g[(2\varepsilon_1-E)+(2\varepsilon_2-E)]=0,
\]

or

\[
E^2-2(\varepsilon_1+\varepsilon_2-g)E
+4\varepsilon_1\varepsilon_2-2g(\varepsilon_1+\varepsilon_2)=0.
\]

This is exactly the determinant of \(H_N^{\rm pair}-EI\), with roots \(E_\pm\). The normalized reciprocal vector has normalization \(\mathcal N^2=\sum_p|2\varepsilon_p-E|^{-2}\). This proof uses only the one-pair sector; it makes no unproved claim about many-pair Bethe completeness or a proposed family of Gaudin charges.

For numerical work, avoid dividing by the small secular denominators. Set \(r=\operatorname{hypot}(d,g)\) and, when \(g>0\),

\[
t=\frac{g/r}{1+|d|/r}.
\]

A lower paired eigenvector is \((1,t)/\sqrt{1+t^2}\) for \(d\le0\), or \((t,1)/\sqrt{1+t^2}\) for \(d>0\). If its coefficients are \((u,v)\), the orthogonal upper vector is \((-v,u)\). Substitute these into the **matrix** eigenvalue equation to verify both branches without using a global coupling variable.

## 4. Singular and degenerate cases

- **\(g=0\), \(d\ne0\):** the paired determinants themselves are eigenvectors. The energies are \(\pm|d|\) for \(H_z\), and \(2\varepsilon_1,2\varepsilon_2\) for \(H_N\). The rational secular equation is inapplicable at these poles; a missing residual is not a successful check.
- **\(d=0\), \(g>0\):** the symmetric paired vector has \(E_z=-2g\), and the antisymmetric vector has \(E_z=0\). The latter has \(c_1+c_2=0\) and \(E_{\rm pair}=2\varepsilon_1\), so the reciprocal Bethe vector is undefined. It remains a valid matrix eigenvector. Together with the four unpaired states, the zero energy has multiplicity five.
- **\(d=g=0\):** \(H_z\) is zero on the entire six-dimensional sector. Any orthonormal basis is valid; comparing one arbitrarily chosen eigenvector with another is not a correctness criterion.
- **Near poles or degeneracy:** raw rational residuals can be badly conditioned even for accurate eigenvectors. Report them honestly, then check a scaled characteristic polynomial, normalization, orthogonality, and matrix residuals. Compare ground eigenspaces rather than particular vectors within a degenerate eigenspace.

## 5. Corrected numerical reference and verification status

For the synthetic parameters \((\varepsilon_1,\varepsilon_2,g)=(-1,-0.5,0.3)\), the Python CLI returns

\[
E_{z,-}=-0.883095189484530,\qquad
E_{{\rm pair},-}=-2.383095189484530.
\]

The previous \(-0.583095189484530\) value omitted the diagonal self-pairing contribution \(-g\). This is a Hamiltonian construction error, not a harmless choice of energy origin: removing the terms changes paired states but leaves the four unpaired states unchanged.

The actual Python run reported energy and full-spectrum errors \(2.220\times10^{-16}\), maximum paired eigenvector residual \(1.570\times10^{-16}\), and passed all ten CLI parameter cases. These are observed floating-point errors for these inputs, not universal bounds at machine epsilon. Reproduce with:

```sh
/usr/bin/python3 -B richardson_h2_verification.py
/usr/bin/python3 -B -m unittest discover -s tests -v
/usr/bin/python3 -B -m experiments.certificates
```

The Python tests independently build CAR operators and compare the complete six-state spectrum, including zero coupling and orbital degeneracy. `experiments/certificate_results.json` separately records synthetic block enclosures, barrier/rate intervals, and fixed-generator CTMC checks, with seeds, tolerances, and units. Dense diagonalization is a finite-size reference, not a scalable chemistry solver.

**Julia status: UNTESTED.** No Julia executable was available in this environment. `richardson_h2_verification.jl` has the same sign, self-pairing, and energy-origin corrections at source level. Its Float64 and optional BigFloat algebra checks must be run separately in a Julia environment; they are not included in the Python success count.

## 6. Limits of the result

A general molecular electronic Hamiltonian contains one- and two-electron terms that this constant-strength pairing model does not represent. Pair correlations being important does not make discarding the other terms exact, including for stretched H₂. To claim a molecular error bound one would need actual integrals, a documented reduction, and a quantitative bound on omitted interactions and other approximations. This derivation proves the stated finite model's spectrum only.
