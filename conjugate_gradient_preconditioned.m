function [x, residuals, iterations] = conjugate_gradient_preconditioned(A, b, x0, M, tol, max_iter)
    % CONJUGATE_GRADIENT_PRECONDITIONED Solve Ax=b using preconditioned Conjugate Gradient
    %
    % Syntax:
    %   [x, residuals, iterations] = conjugate_gradient_preconditioned(A, b, x0, M, tol, max_iter)
    %
    % Input:
    %   A         - Coefficient matrix (symmetric positive definite)
    %   b         - Right-hand side vector
    %   x0        - Initial solution guess
    %   M         - Preconditioner matrix (M = L*L' from incomplete Cholesky)
    %   tol       - Convergence tolerance for residual norm
    %   max_iter  - Maximum number of iterations
    %
    % Output:
    %   x         - Solution vector
    %   residuals - Norm of residual at each iteration
    %   iterations- Number of iterations performed
    %
    % The preconditioner M is applied as: M^{-1} is implicitly computed
    % using the factorization M = L*L' and solving L*y = r via back substitution

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

    % Apply preconditioner: y0 = M^{-1} * r0
    % Since M = L*L', solve L*z = r and then L'*y = z
    z = M \ r;  % Solve L*z = r (back substitution from ichol)
    y = M' \ z; % Solve L'*y = z

    % Initial search direction: p0 = y0
    p = y;

    rTy = r' * y;  % Store initial r'*y for beta computation

    for k = 1:max_iter
        % Compute step size: alpha_k = (r_k^T * y_k) / (p_k^T * A * p_k)
        Ap = A * p;
        pTAp = p' * Ap;

        alpha = rTy / pTAp;

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

        % Apply preconditioner: y_{k+1} = M^{-1} * r_{k+1}
        z = M \ r;
        y_new = M' \ z;

        % Compute beta_k = (r_{k+1}^T * y_{k+1}) / (r_k^T * y_k)
        rTy_new = r' * y_new;
        beta = rTy_new / rTy;

        % Update search direction: p_{k+1} = y_{k+1} + beta_k * p_k
        p = y_new + beta * p;

        % Store r'*y for next iteration
        y = y_new;
        rTy = rTy_new;
    end

    iterations = max_iter;
end
