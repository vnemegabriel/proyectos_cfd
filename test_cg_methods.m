%% Unit Test: Conjugate Gradient Methods
% Verify correctness of CG implementation and preconditioned variants

clear all; close all; clc;

fprintf('\n========== UNIT TESTS: Conjugate Gradient Methods ==========\n\n');

%% Test 1: Small system verification
fprintf('Test 1: Small system (3x3 connectivity matrix)...\n');

A_test = generate_connectivity_matrix(3);
n_test = 9;

% Create right-hand side
b_test = ones(n_test, 1);
x0_test = zeros(n_test, 1);

% Solve with CG
tol_test = 1e-10;
[x_cg_test, res_test, iter_test] = conjugate_gradient_method(A_test, b_test, x0_test, tol_test, n_test);

% Verify solution
residual = A_test * x_cg_test - b_test;
error = norm(residual);

fprintf('  - System size: %d x %d\n', n_test, n_test);
fprintf('  - Iterations to convergence: %d\n', iter_test);
fprintf('  - Final residual norm: %.6e\n', norm(res_test(end)));
fprintf('  - Solution error ||Ax - b||: %.6e\n', error);

if error < tol_test * 10
    fprintf('  ✓ PASSED: Solution is accurate\n\n');
else
    fprintf('  ✗ FAILED: Solution error too large\n\n');
end

%% Test 2: Symmetric positive definiteness
fprintf('Test 2: Verify SPD property of connectivity matrix...\n');

A = generate_connectivity_matrix(5);
eigs = eig(full(A));
is_symmetric = isequal(A, A');
all_positive = all(eigs > 0);

fprintf('  - Matrix symmetric: %s\n', string(is_symmetric));
fprintf('  - All eigenvalues > 0: %s\n', string(all_positive));
fprintf('  - Min eigenvalue: %.6e\n', min(eigs));
fprintf('  - Max eigenvalue: %.6e\n', max(eigs));

if is_symmetric && all_positive
    fprintf('  ✓ PASSED: Matrix is SPD\n\n');
else
    fprintf('  ✗ FAILED: Matrix is not SPD\n\n');
end

%% Test 3: Preconditioner quality
fprintf('Test 3: Incomplete Cholesky factorization quality...\n');

A = generate_connectivity_matrix(5);
L = ichol(A);

% Reconstruct M = L*L'
M = L * L';

% Compare with original matrix
rel_error = norm(M - A, 'fro') / norm(A, 'fro');
nnz_A = nnz(A);
nnz_L = nnz(L);

fprintf('  - Original matrix A: %d non-zeros\n', nnz_A);
fprintf('  - Incomplete Cholesky L: %d non-zeros\n', nnz_L);
fprintf('  - Relative error ||L*L'' - A||_F / ||A||_F: %.6e\n', rel_error);

if rel_error < 0.1  % Incomplete Cholesky won't be exact
    fprintf('  ✓ PASSED: Good approximation\n\n');
else
    fprintf('  ⚠ WARNING: Approximation quality may be limited\n\n');
end

%% Test 4: CG vs Preconditioned CG on small system
fprintf('Test 4: Comparison of CG variants on m=5 system...\n');

A = generate_connectivity_matrix(5);
b = (1/25) * ones(25, 1);
x0 = ones(25, 1);
tol = 1e-12;
max_iter = 25;

% Solve with CG
[x_cg, res_cg, iter_cg] = conjugate_gradient_method(A, b, x0, tol, max_iter);

% Solve with preconditioned CG
L = ichol(A);
[x_pcg, res_pcg, iter_pcg] = conjugate_gradient_preconditioned(A, b, x0, L, tol, max_iter);

fprintf('  - Unpreconditioned CG iterations: %d\n', iter_cg);
fprintf('  - Preconditioned CG iterations: %d\n', iter_pcg);
fprintf('  - Improvement factor: %.2f×\n', iter_cg/iter_pcg);

% Both should reach similar solution accuracy
x_exact = A \ b;
error_cg = norm(x_cg - x_exact) / norm(x_exact);
error_pcg = norm(x_pcg - x_exact) / norm(x_exact);

fprintf('  - CG relative error: %.6e\n', error_cg);
fprintf('  - PCG relative error: %.6e\n', error_pcg);

if iter_pcg < iter_cg && error_pcg < 1e-10 && error_cg < 1e-10
    fprintf('  ✓ PASSED: Preconditioner improves convergence\n\n');
else
    fprintf('  ✗ FAILED: Unexpected behavior\n\n');
end

%% Test 5: Convergence tolerance
fprintf('Test 5: Convergence tolerance behavior...\n');

A = generate_connectivity_matrix(4);
b = 0.0625 * ones(16, 1);
x0 = ones(16, 1);

tolerances = [1e-4, 1e-6, 1e-8, 1e-10, 1e-12];
iters = [];

for tol = tolerances
    [~, ~, iter] = conjugate_gradient_method(A, b, x0, tol, 16);
    iters = [iters; iter];
end

fprintf('  Tolerance  | Iterations\n');
fprintf('  -----------+----------\n');
for i = 1:length(tolerances)
    fprintf('  %.6e |  %d\n', tolerances(i), iters(i));
end

if all(diff(iters) >= 0)
    fprintf('  ✓ PASSED: Iterations increase with stricter tolerance\n\n');
else
    fprintf('  ✗ FAILED: Unexpected iteration behavior\n\n');
end

fprintf('========== ALL TESTS COMPLETED ==========\n\n');
