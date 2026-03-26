%% Exercise 7: Conjugate Gradient Method with Preconditioners
% Computational Matrix Algebra - Dr. Santiago Márquez Damián
%
% This script solves Exercise 7 from AlgebraMatComp:
% "Implement the Conjugate Gradient Method and extend it to use preconditioners
%  (Incomplete Cholesky and Modified Incomplete Cholesky). Compare residual
%  evolution and draw conclusions."

clear all; close all; clc;

%% Setup Parameters
fprintf('\n========== Exercise 7: CG with Preconditioners ==========\n\n');

% Matrix size parameter: m x m grid
m = 20;
n = m * m;

% Create connectivity matrix (5-point Laplacian stencil)
fprintf('Creating %dx%d connectivity matrix...\n', m, m);
A = generate_connectivity_matrix(m);

% Verify that A is symmetric positive definite (SPD)
eigenvalues = eig(full(A));
is_symmetric = isequal(A, A');
is_spd = all(eigenvalues > 0);

fprintf('Matrix Properties:\n');
fprintf('  - Size: %d x %d (%d unknowns)\n', n, n, n);
fprintf('  - Symmetric: %s\n', string(is_symmetric));
fprintf('  - Positive Definite: %s\n', string(is_spd));
fprintf('  - Min eigenvalue: %.6e\n', min(eigenvalues));
fprintf('  - Max eigenvalue: %.6e\n', max(eigenvalues));
fprintf('  - Condition number (κ): %.6e\n\n', max(eigenvalues) / min(eigenvalues));

%% Setup System of Equations: Ax = b
% Parameters: scale factor = 1/(m*m), b = 1/(m*m), x0 = 1
scale_factor = 1 / (m * m);
b = scale_factor * ones(n, 1);
x0 = ones(n, 1);

% Solve for exact solution using MATLAB backslash
x_exact = A \ b;

fprintf('System Setup:\n');
fprintf('  - Scale factor: %.6e\n', scale_factor);
fprintf('  - Right-hand side b: all %.6e\n', scale_factor);
fprintf('  - Initial guess x0: all 1.0\n');
fprintf('  - Convergence tolerance: 1e-12\n\n');

%% Solver Parameters
tol = 1e-12;
max_iter = n;  % Maximum iterations

%% Solver 1: Conjugate Gradient (no preconditioner)
fprintf('Solver 1: CG without preconditioner...\n');
tic;
[x_cg, res_cg, iter_cg] = conjugate_gradient_method(A, b, x0, tol, max_iter);
time_cg = toc;

error_cg = norm(x_cg - x_exact) / norm(x_exact);
fprintf('  - Iterations: %d\n', iter_cg);
fprintf('  - Time: %.6f seconds\n', time_cg);
fprintf('  - Final residual: %.6e\n', res_cg(end));
fprintf('  - Relative error: %.6e\n\n', error_cg);

%% Solver 2: Preconditioned CG with Incomplete Cholesky
fprintf('Solver 2: CG with Incomplete Cholesky preconditioner...\n');
tic;
% Generate incomplete Cholesky factorization: A ≈ L*L'
% ichol returns the lower triangular factor L
L_ichol = ichol(A);
time_setup_ic = toc;

tic;
[x_ic, res_ic, iter_ic] = conjugate_gradient_preconditioned(A, b, x0, L_ichol, tol, max_iter);
time_ic = toc;

error_ic = norm(x_ic - x_exact) / norm(x_exact);
fprintf('  - Setup time (ichol): %.6f seconds\n', time_setup_ic);
fprintf('  - Solver time: %.6f seconds\n', time_ic);
fprintf('  - Total time: %.6f seconds\n', time_setup_ic + time_ic);
fprintf('  - Iterations: %d\n', iter_ic);
fprintf('  - Final residual: %.6e\n', res_ic(end));
fprintf('  - Relative error: %.6e\n\n', error_ic);

%% Solver 3: Preconditioned CG with Modified Incomplete Cholesky
fprintf('Solver 3: CG with Modified Incomplete Cholesky preconditioner...\n');
tic;
% Modified incomplete Cholesky with relaxation parameter (default 0.95)
options.type = 'nofill';
options.droptol = 0;  % No dropping
options.michol = 'on'; % Modified incomplete Cholesky
L_michol = ichol(A, options);
time_setup_mic = toc;

tic;
[x_mic, res_mic, iter_mic] = conjugate_gradient_preconditioned(A, b, x0, L_michol, tol, max_iter);
time_mic = toc;

error_mic = norm(x_mic - x_exact) / norm(x_exact);
fprintf('  - Setup time (michol): %.6f seconds\n', time_setup_mic);
fprintf('  - Solver time: %.6f seconds\n', time_mic);
fprintf('  - Total time: %.6f seconds\n', time_setup_mic + time_mic);
fprintf('  - Iterations: %d\n', iter_mic);
fprintf('  - Final residual: %.6e\n', res_mic(end));
fprintf('  - Relative error: %.6e\n\n', error_mic);

%% Comparison Summary
fprintf('========== COMPARISON SUMMARY ==========\n\n');
fprintf('%-30s | %8s | %8s | %12s\n', 'Method', 'Iterations', 'Time (s)', 'Error');
fprintf('%-30s | %8s | %8s | %12s\n', repmat('-', 1, 30), repmat('-', 1, 8), repmat('-', 1, 8), repmat('-', 1, 12));
fprintf('%-30s | %8d | %8.6f | %12.6e\n', 'CG (no precond.)', iter_cg, time_cg, error_cg);
fprintf('%-30s | %8d | %8.6f | %12.6e\n', 'CG + Incomplete Cholesky', iter_ic, time_setup_ic + time_ic, error_ic);
fprintf('%-30s | %8d | %8.6f | %12.6e\n', 'CG + Mod. Incomplete Cholesky', iter_mic, time_setup_mic + time_mic, error_mic);

% Compute speedup relative to CG
speedup_ic = iter_cg / iter_ic;
speedup_mic = iter_cg / iter_mic;
fprintf('\nIteration Reduction:\n');
fprintf('  - Incomplete Cholesky: %.2f× fewer iterations\n', speedup_ic);
fprintf('  - Modified Inc. Cholesky: %.2f× fewer iterations\n', speedup_mic);

%% Visualization: Residual Evolution
figure('Position', [100 100 1200 500]);

% Plot 1: Linear scale
subplot(1, 2, 1);
iterations_cg = 0:length(res_cg)-1;
iterations_ic = 0:length(res_ic)-1;
iterations_mic = 0:length(res_mic)-1;

semilogy(iterations_cg, res_cg, 'b-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'CG (no precond.)');
hold on;
semilogy(iterations_ic, res_ic, 'r-s', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'CG + Incomplete Cholesky');
semilogy(iterations_mic, res_mic, 'g-^', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'CG + Mod. Inc. Cholesky');

grid on;
xlabel('Iteration', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('||r||_2', 'FontSize', 12, 'FontWeight', 'bold');
title('Residual Evolution (Linear Iteration Scale)', 'FontSize', 13, 'FontWeight', 'bold');
legend('FontSize', 10, 'Location', 'best');
set(gca, 'YScale', 'log');

% Plot 2: Logarithmic scale (log-log to see convergence rate)
subplot(1, 2, 2);
iterations_cg_plot = iterations_cg(iterations_cg > 0);
iterations_ic_plot = iterations_ic(iterations_ic > 0);
iterations_mic_plot = iterations_mic(iterations_mic > 0);

loglog(iterations_cg_plot, res_cg(iterations_cg_plot+1), 'b-o', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'CG (no precond.)');
hold on;
loglog(iterations_ic_plot, res_ic(iterations_ic_plot+1), 'r-s', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'CG + Incomplete Cholesky');
loglog(iterations_mic_plot, res_mic(iterations_mic_plot+1), 'g-^', 'LineWidth', 2, 'MarkerSize', 4, 'DisplayName', 'CG + Mod. Inc. Cholesky');

grid on;
xlabel('Iteration (log scale)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('||r||_2 (log scale)', 'FontSize', 12, 'FontWeight', 'bold');
title('Residual Evolution (Log-Log Scale)', 'FontSize', 13, 'FontWeight', 'bold');
legend('FontSize', 10, 'Location', 'best');

sgtitle(sprintf('Conjugate Gradient Method with Preconditioners (m=%d, n=%d)', m, n), ...
    'FontSize', 14, 'FontWeight', 'bold');

savefig('exercise7_residual_evolution.fig');
saveas(gcf, 'exercise7_residual_evolution.png');
fprintf('\n✓ Figure saved as "exercise7_residual_evolution.fig" and "exercise7_residual_evolution.png"\n');

%% Conclusions
fprintf('\n========== CONCLUSIONS ==========\n\n');
fprintf('1. CONVERGENCE IMPROVEMENT:\n');
fprintf('   Both incomplete Cholesky variants significantly reduce the number of\n');
fprintf('   iterations required for convergence compared to the unpreconditioned CG.\n');
fprintf('   - Unpreconditioned CG: %d iterations\n', iter_cg);
fprintf('   - With Incomplete Cholesky: %d iterations (%.1f%% reduction)\n', iter_ic, (1-iter_ic/iter_cg)*100);
fprintf('   - With Mod. Inc. Cholesky: %d iterations (%.1f%% reduction)\n\n', iter_mic, (1-iter_mic/iter_cg)*100);

fprintf('2. PRECONDITIONER EFFECTIVENESS:\n');
fprintf('   The preconditioner works by approximating the inverse of A,\n');
fprintf('   improving the spectral properties of the system and reducing\n');
fprintf('   the effective condition number of the preconditioned system.\n\n');

fprintf('3. MODIFIED VS. INCOMPLETE CHOLESKY:\n');
if iter_mic <= iter_ic
    fprintf('   Modified Incomplete Cholesky requires fewer or equal iterations\n');
    fprintf('   due to better approximation of the original matrix factorization.\n');
else
    fprintf('   Incomplete Cholesky achieves comparable results with simpler setup.\n');
end
fprintf('\n');

fprintf('4. PRACTICAL IMPLICATIONS:\n');
fprintf('   - Setup cost: One-time factorization (ichol)\n');
fprintf('   - Solver cost: Each preconditioned iteration requires solving\n');
fprintf('     L*y = r and L^T*z = y (via back-substitution, O(nnz(L)) per iteration)\n');
fprintf('   - Net benefit: Significant reduction in total computational cost\n');
fprintf('   - Trade-off: Setup time vs. iteration reduction\n\n');

fprintf('5. MATRIX CHARACTERISTICS:\n');
fprintf('   This sparse, symmetric positive definite matrix appears in:\n');
fprintf('   - Finite Element Method (FEM) for elliptic PDEs\n');
fprintf('   - Finite Difference Method (FDM) with uniform grids\n');
fprintf('   - Graph Laplacians in signal processing and machine learning\n');
fprintf('   - The 5-point stencil is the discrete Laplacian operator in 2D\n\n');

fprintf('========================================\n\n');
