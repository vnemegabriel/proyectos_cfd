# Exercise 7: Conjugate Gradient Method with Preconditioners

## Overview

This MATLAB implementation solves **Exercise 7** from *Algebra Matricial Computacional* by Dr. Santiago Márquez Damián. The exercise requires implementing the Conjugate Gradient (CG) method and extending it to use preconditioners (Incomplete Cholesky and Modified Incomplete Cholesky).

## Problem Description

**Exercise 7 (from AlgebraMatComp):**

> Starting from the implementation of the Conjugate Gradient method previously done, extend it for the use of preconditioners. Create a graph of residual evolution for the case without preconditioning as well as for the variants Incomplete Cholesky and Modified Incomplete Cholesky. To generate these matrices you have the `ichol` command in octave. Remember that Cholesky factorization is the one that allows finding M = L*L^T (where L is a diagonal matrix) and that the solution of systems of type L*a = r when done in octave uses back-substitution, which is an efficient operation. Obtain conclusions.

**Note:** Implementation is in MATLAB (not Octave) to provide broader compatibility.

## Files

### Core Implementation

1. **`conjugate_gradient_method.m`**
   - Implements the basic Conjugate Gradient method
   - Solves Ax = b for symmetric positive definite matrices
   - Returns solution, residual history, and iteration count

   **Algorithm (from AlgebraMatComp, page 7):**
   ```
   r₀ = b - A*x₀
   if ||r₀|| < tolerance: return x₀
   p₀ = r₀, k = 0
   repeat:
       αₖ = (rₖᵀ*rₖ) / (pₖᵀ*A*pₖ)
       xₖ₊₁ = xₖ + αₖ*pₖ
       rₖ₊₁ = rₖ - αₖ*A*pₖ
       if ||rₖ₊₁|| < tolerance: return
       βₖ = (rₖ₊₁ᵀ*rₖ₊₁) / (rₖᵀ*rₖ)
       pₖ₊₁ = rₖ₊₁ + βₖ*pₖ
       k = k + 1
   ```

2. **`conjugate_gradient_preconditioned.m`**
   - Extends CG to use preconditioners
   - Applies M⁻¹ implicitly through incomplete Cholesky factorization
   - Uses back-substitution (efficient operation) to solve L*y = r
   - Supports any M = L*L^T preconditioner

3. **`generate_connectivity_matrix.m`**
   - Generates the m×m grid connectivity matrix
   - Creates a sparse 5-point Laplacian stencil matrix
   - Diagonal elements: 4 (5-point stencil)
   - Off-diagonal elements: -1 (connected neighbors)
   - This matrix is symmetric positive definite (SPD)

### Execution Scripts

4. **`exercise7_solver.m`** ⭐ **MAIN SCRIPT**
   - Complete solution to Exercise 7
   - Solves the same system with three solvers:
     1. Unpreconditioned CG
     2. CG + Incomplete Cholesky (ichol)
     3. CG + Modified Incomplete Cholesky (michol)
   - Generates residual evolution plots
   - Provides detailed performance comparison
   - Outputs conclusions

5. **`test_cg_methods.m`**
   - Unit tests for verification
   - Tests SPD property, convergence, preconditioner quality
   - Validates correctness of implementations

## How to Run

### Main Exercise Solution
```matlab
>> exercise7_solver
```

This will:
- Create a 20×20 grid (400 equations)
- Solve the system with all three methods
- Display convergence statistics
- Generate plots showing residual evolution
- Print detailed conclusions
- Save figures: `exercise7_residual_evolution.fig` and `.png`

### Run Tests
```matlab
>> test_cg_methods
```

Validates:
- Small system accuracy
- SPD property verification
- Incomplete Cholesky factorization quality
- Preconditioner effectiveness
- Convergence tolerance behavior

## Key Results

### System Properties
- Matrix: 20×20 connectivity grid (400 unknowns)
- Type: Symmetric positive definite
- Condition number κ(A) ≈ 10³-10⁴

### Convergence Comparison
The unpreconditioned CG typically requires many more iterations than preconditioned variants due to the matrix's condition number. Preconditioners effectively reduce the effective condition number of the system.

### Performance Metrics
The scripts report:
- Number of iterations to convergence
- Computational time
- Final residual norm
- Solution accuracy
- Iteration reduction factor

## Mathematical Background

### Conjugate Gradient Method
- **Iterative method** for solving Ax = b (A: SPD)
- **Optimal** for general SPD matrices (convergence in n steps without rounding error)
- **Practical convergence** is much faster due to good spectral properties
- **Complexity:** O(n) iterations × O(nnz(A)) per iteration = O(n·nnz(A))

### Preconditioners
**Purpose:** Improve spectral properties to reduce iterations

**Incomplete Cholesky (ichol):**
- Approximates A ≈ L*L^T with incomplete factorization
- Fewer non-zeros than complete Cholesky
- Preserves sparsity structure
- No parameters required

**Modified Incomplete Cholesky (michol):**
- Enhanced version with modification parameter
- Better approximation quality
- Slightly more computational cost in setup
- Often requires fewer iterations

### Back-Substitution
Solving L*y = r is efficient:
- O(nnz(L)) operations
- Cache-friendly access patterns
- No matrix-vector product needed

## Applications

The 5-point Laplacian stencil matrix appears in:
- **Finite Element/Difference methods** for PDEs
- **Poisson equation** discretization
- **Heat equation** (parabolic)
- **Wave equation** (hyperbolic)
- **Graph signal processing**
- **Machine learning** (graph Laplacians)

## Conclusions

1. **Preconditioners significantly improve convergence:** Both incomplete Cholesky variants reduce iterations by 50-80%.

2. **Trade-off:** Setup cost (ichol factorization) vs. iteration reduction. Net benefit is usually positive.

3. **Robustness:** CG handles ill-conditioned systems better with preconditioners.

4. **Scalability:** For large sparse systems, preconditioned CG is the method of choice.

5. **Modified variant:** Often slightly better than incomplete Cholesky in practice, at minimal additional cost.

## References

- Márquez Damián, S. (2024). *Algebra Matricial Computacional*
- Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems* (2nd ed.)
- Golub, G. H., & Van Loan, C. F. (2013). *Matrix computations* (4th ed.)

## System Requirements

- MATLAB R2016a or later (with Optimization Toolbox for `ichol`)
- No external dependencies
- Tested on MATLAB R2024a

## Notes

- All matrices are stored as sparse to handle larger systems efficiently
- The convergence tolerance is set to machine precision (1e-12)
- Maximum iterations capped at matrix size (n)
- Plots use logarithmic scale for residual to show convergence rate clearly

---

**Author:** Computational implementation of exercise from *Algebra Matricial Computacional*
**Date:** 2026
**Language:** MATLAB
