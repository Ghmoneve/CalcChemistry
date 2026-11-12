# CalcChemistry: Formal Verification of Quantum Chemistry Theory

## Goal

Formalize mathematical foundations for polynomial-time molecular calculations using Lean 4.

## Research Focus

Explore whether specific molecular Hamiltonian structures admit exact solutions that bypass exponential complexity:
- Quantum integrable systems (Bethe ansatz, Yang-Baxter)
- Reduced density matrix methods with tight algebraic-geometric bounds
- Tensor network decompositions with optimal rank bounds
- Lie algebra and quantum group symmetries
- Rigorous complexity analysis of physical approximations

## Structure

- `CalcChemistry/Basic.lean` - Core definitions (Hilbert spaces, Hamiltonians, Born-Oppenheimer)
- `CalcChemistry/Integrable.lean` - (planned) Integrable system characterization
- `CalcChemistry/RDM.lean` - (planned) N-representability and SDP hierarchy
- `CalcChemistry/TensorNetwork.lean` - (planned) Tensor decomposition theorems
- `CalcChemistry/Complexity.lean` - (planned) Formal complexity bounds

## Build

```bash
lake build
```

## Citation Policy

All theorems are original or cite sources explicitly. No fabricated references.
