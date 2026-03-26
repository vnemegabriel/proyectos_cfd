function [x, residuals, iterations] = conjugate_gradient_method(A, b, x0, tol, max_iter)
    % CONJUGATE_GRADIENT_METHOD Solve Ax = b using Conjugate Gradient Method
    %
    % Syntax:
    %   [x, residuals, iterations] = conjugate_gradient_method(A, b, x0, tol, max_iter)
    %
    % Input:
    %   A         - Coefficient matrix (symmetric positive definite)
    %   b         - Right-hand side vector
    %   x0        - Initial solution guess
    %   tol       - Convergence tolerance for residual norm
    %   max_iter  - Maximum number of iterations
    %
    % Output:
    %   x         - Solution vector
    %   residuals - Norm of residual at each iteration
    %   iterations- Number of iterations performed
    %
    % Reference: Algorithm from "Algebra matricial computacional" by Dr. Santiago Márquez Damián

    n = length(b);
    x = x0;

    % Initial residual: r0 = b - A*x0
    r = b - A * x;
    r_norm = norm(r);

    residuals = [r_norm];

    % Check if initial solution is already good enough
    if r_norm < tol
        iterations = 0;
        return;
    end

    % Initial search direction: p0 = r0
    p = r;

    for k = 1:max_iter
        % Compute step size alpha_k = (r_k^T * r_k) / (p_k^T * A * p_k)
        rTr = r' * r;
        Ap = A * p;
        pTAp = p' * Ap;

        alpha = rTr / pTAp;

        % Update solution: x_{k+1} = x_k + alpha_k * p_k
        x = x + alpha * p;

        % Update residual: r_{k+1} = r_k - alpha_k * A * p_k
        r = r - alpha * Ap;
        r_norm = norm(r);

        residuals = [residuals; r_norm];

        % Check convergence
        if r_norm < tol
            iterations = k;
            return;
        end

        % Compute beta_k = (r_{k+1}^T * r_{k+1}) / (r_k^T * r_k)
        rTr_new = r' * r;
        beta = rTr_new / rTr;

        % Update search direction: p_{k+1} = r_{k+1} + beta_k * p_k
        p = r + beta * p;
    end

    iterations = max_iter;
end
