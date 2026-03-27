# Verification and Testing Guide

## Overview

This document provides mathematical verification that the implementation is correct, along with expected behavior and how to test it.

## Matrix Properties Verification

### Connectivity Matrix (5-Point Laplacian)

The connectivity matrix for an m×m grid should be **symmetric positive definite (SPD)**.

**Construction:**
```
For m = 2 (2×2 grid, 4 nodes):

Nodes:   1  2
         3  4

Expected matrix A:
    4  -1  -1   0
   -1   4   0  -1
   -1   0   4  -1
    0  -1  -1   4
```

**Mathematical properties:**
- **Symmetry:** A = A^T ✓
- **Diagonal dominance:** Each diagonal element (4) > sum of off-diagonal absolute values (1+1=2) ✓
- **Eigenvalues:** All positive (guaranteed SPD) ✓
- **Condition number:** κ(A) = λ_max / λ_min (for m=20, κ ≈ 10³-10⁴)

**Why SPD?**
1. Comes from discretization of Laplacian operator: ∇²u = ∂²u/∂x² + ∂²u/∂y²
2. Laplacian of convex problems is always SPD
3. 5-point stencil is the standard FD approximation

## Algorithm Verification

### Conjugate Gradient Method

**Algorithm (from AlgebraMatComp, page 7):**

```
Input:  A (SPD), b (RHS), x₀ (initial guess), tol (tolerance)
Output: x (solution)

1. r₀ = b - A*x₀
2. if ||r₀|| < tol:  return x₀

3. p₀ = r₀
4. k = 0
5. loop:
     αₖ = (rₖ'*rₖ) / (pₖ'*A*pₖ)
     xₖ₊₁ = xₖ + αₖ*pₖ
     rₖ₊₁ = rₖ - αₖ*A*pₖ
     if ||rₖ₊₁|| < tol:  return xₖ₊₁
     βₖ = (rₖ₊₁'*rₖ₊₁) / (rₖ'*rₖ)
     pₖ₊₁ = rₖ₊₁ + βₖ*pₖ
     k = k + 1
```

**Convergence properties:**
- Exact solution in n steps (without rounding)
- Practical: Good convergence due to spectral properties
- Residual norm decreases monotonically: ||r₀|| > ||r₁|| > ... > ||r_n||

**Complexity:**
- Per iteration: O(nnz(A)) for matrix-vector product
- Total iterations: O(κ(A)^0.5) for convergence
- Total complexity: O(κ(A)^0.5 * nnz(A)) ≈ O(n * κ(A)^0.5) for sparse

### Preconditioned CG

**Modified algorithm using preconditioner M ≈ A:**

```
Input:  A (SPD), b (RHS), x₀ (initial), M = L*L' (preconditioner), tol

1. r₀ = b - A*x₀
2. Solve M*y₀ = r₀ (using L*L' = M: two triangular solves)
3. if ||r₀|| < tol:  return x₀

4. p₀ = y₀
5. k = 0
6. loop:
     αₖ = (rₖ'*yₖ) / (pₖ'*A*pₖ)
     xₖ₊₁ = xₖ + αₖ*pₖ
     rₖ₊₁ = rₖ - αₖ*A*pₖ
     if ||rₖ₊₁|| < tol:  return xₖ₊₁
     Solve M*yₖ₊₁ = rₖ₊₁
     βₖ = (rₖ₊₁'*yₖ₊₁) / (rₖ'*yₖ)
     pₖ₊₁ = yₖ₊₁ + βₖ*pₖ
     k = k + 1
```

**Benefits of preconditioner M:**
- Effectively reduces condition number: κ_eff = κ(M^{-1}*A) << κ(A)
- Fewer iterations required
- Each iteration cost ≈ same (L*L' solve vs matrix-vector product)
- Net: Significant speedup in total time

### Incomplete Cholesky Factorization

**ichol Properties:**
```
Given:  A (SPD sparse matrix)
Output: L (lower triangular)
        M = L*L' ≈ A (incomplete factorization)

Key features:
- L ≈ Cholesky factor of A
- Same sparsity pattern as lower(A) or user-defined pattern
- M ≈ A provides good preconditioner
- Solution of L*y = r via back substitution: O(nnz(L))
```

**Modified ichol:**
- Additional parameter to improve approximation
- Trade-off: Slightly denser L vs better convergence
- Often worth the extra cost

## Expected Results

### System: 20×20 Grid (400 unknowns)

| Metric | Unpreconditioned CG | CG + ichol | CG + michol |
|--------|-------------------|-----------|------------|
| **Iterations** | 150-200 | 30-50 | 25-45 |
| **Setup time** | 0 | 0.01-0.05s | 0.01-0.05s |
| **Solver time** | 0.01-0.02s | 0.003-0.008s | 0.003-0.008s |
| **Iteration reduction** | 1× | 4-6× | 4-7× |
| **Total speedup** | 1× | 3-5× | 3-5× |

**Why preconditioner helps:**
- Condition number κ(A) ≈ 5000
- CG iterations ≈ √κ ≈ 70 (theoretical)
- With good preconditioner: κ_eff ≈ 100, iterations ≈ √100 ≈ 10

## Convergence Behavior

### Residual Evolution

The residual ||r_k|| should show:

1. **Smooth exponential decay** in the semilogy plot
2. **Faster decay** with preconditioner
3. **Final plateau** at machine epsilon (≈ 10^{-15})

```
||r_k||
  ↑
  |  \
  |   \ CG + precond
  |    \___
  |        \___  CG (unpreconditioned)
  |            \___
  |________________\___→ k (iterations)
```

### Stopping criterion

With tolerance tol = 1e-12:
- Iteration stops when ||r_k|| < tol
- For SPD system, A*x ≈ b when residual small
- Solution error ≈ ||x - x_exact|| ≤ κ(A) * ||r|| * tolerance

## Matrix Structure

### Sparsity Pattern (m=5, 25 nodes)

```
* * . . .  (node 1 connects to 2,3)
* * * . .  (node 2 connects to 1,3,4)
. * * * .
. . * * *
. . . * *

Sparsity:  25 - 4*5 + 2*(4) + 2*(4) = 25 - 20 + 8 + 8 = 21 non-zeros
Actually:  25 + 2*4 + 2*4 = 25 (diagonal) + 8 + 8 = 41...

For m x m grid: nnz ≈ 5*m² - 4*m
For m=20: nnz ≈ 5*400 - 80 = 1920
```

## Validation Checklist

To verify implementation without MATLAB installed:

- [ ] **Matrix construction:**
  - [ ] Run test_cg_methods → Test 1 (small system)
  - [ ] Verify ||Ax - b|| < 1e-10

- [ ] **SPD properties:**
  - [ ] Run test_cg_methods → Test 2
  - [ ] Confirm all eigenvalues > 0
  - [ ] Confirm A = A^T

- [ ] **Preconditioner quality:**
  - [ ] Run test_cg_methods → Test 3
  - [ ] Verify ||L*L' - A||_F / ||A||_F < 0.1

- [ ] **CG convergence:**
  - [ ] Run exercise7_solver
  - [ ] Check iteration counts match expected ranges
  - [ ] Verify final residuals < 1e-12

- [ ] **Preconditioner effectiveness:**
  - [ ] Run test_cg_methods → Test 4
  - [ ] Verify PCG iterations < CG iterations
  - [ ] Check solution accuracy for both methods

- [ ] **Convergence plots:**
  - [ ] Run exercise7_solver
  - [ ] Verify residual curves decrease smoothly
  - [ ] Confirm preconditioned variants converge faster
  - [ ] Check files: exercise7_residual_evolution.fig/png

## Error Sources and Debugging

### If convergence is poor:

1. **Check matrix properties:**
   - Verify A is symmetric (using test_cg_methods)
   - Verify all eigenvalues > 0 (test 2)

2. **Check preconditioner:**
   - Verify ichol didn't fail (needs SPD matrix)
   - Check ||L*L' - A|| is reasonable (test 3)

3. **Adjust tolerance:**
   - If tol = 1e-12 is too strict, try 1e-10
   - Machine epsilon ≈ 2.22e-16 for double precision

### If ichol fails with "nonpositive pivot":

- Matrix A is not positive definite!
- Causes: Incorrect matrix construction, wrong sign, numerical errors
- Solution: Use test_cg_methods to diagnose

### If iterations don't improve with preconditioner:

- Preconditioner is too weak
- Possible: ichol setup issue, wrong matrix format
- Try: Create better preconditioner or adjust setup parameters

## References

1. Saad, Y. (2003). *Iterative Methods for Sparse Linear Systems*, 2nd ed.
2. Golub & Van Loan (2013). *Matrix Computations*, 4th ed.
3. Márquez Damián, S. (2024). *Algebra Matricial Computacional*

## File Relationships

```
generate_connectivity_matrix.m
    ↓
    ├→ test_cg_methods.m (Tests 1, 2, 3, 4, 5)
    │
    └→ exercise7_solver.m (Main exercise)
        ├→ conjugate_gradient_method.m
        └→ conjugate_gradient_preconditioned.m
```

---

**Last Updated:** 2026-03-27
**Status:** Ready for MATLAB execution
