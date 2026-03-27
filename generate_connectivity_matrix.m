function A = generate_connectivity_matrix(m)
    % GENERATE_CONNECTIVITY_MATRIX Create connectivity matrix for m x m grid
    %
    % Syntax:
    %   A = generate_connectivity_matrix(m)
    %
    % Input:
    %   m - Size of grid (m x m)
    %
    % Output:
    %   A - Connectivity matrix (n x n where n = m*m)
    %       Symmetric positive definite 5-point Laplacian stencil
    %       - Diagonal: 4 (number of neighbors in 5-point stencil)
    %       - Off-diagonal: -1 (for adjacent grid cells)
    %
    % Description:
    %   Creates a sparse connectivity matrix for m x m grid discretization.
    %   This is the classic 5-point Laplacian stencil used in:
    %   - Finite Element Method (FEM) for Poisson problems
    %   - Finite Difference Method (FDM) on rectangular grids
    %   - Graph Laplacian matrices
    %
    % Reference: Exercise 2 from "Algebra matricial computacional"
    %            by Dr. Santiago Márquez Damián
    %
    % Theory: The 5-point stencil approximates the 2D Laplacian:
    %   -∆u ≈ 4*u_c - u_north - u_south - u_east - u_west
    %
    % The resulting matrix is symmetric positive definite (SPD).

    n = m * m;

    % Pre-allocate storage for sparse matrix entries
    max_entries = 3 * n;  % n diagonal + 2*(n-m) horizontal + 2*(n-m) vertical
    rows = zeros(max_entries, 1);
    cols = zeros(max_entries, 1);
    vals = zeros(max_entries, 1);
    ent = 0;  % Entry counter

    % Build connectivity matrix using 5-point stencil (lower triangle only)
    for i = 1:m
        for j = 1:m
            % Node index in linear ordering (by rows)
            node = (i - 1) * m + j;

            % Diagonal element: 4 (self-connectivity strength)
            ent = ent + 1;
            rows(ent) = node;
            cols(ent) = node;
            vals(ent) = 4;

            % Right neighbor (same row, next column) - add to lower triangle
            if j < m
                neighbor = (i - 1) * m + (j + 1);
                ent = ent + 1;
                rows(ent) = neighbor;
                cols(ent) = node;
                vals(ent) = -1;
            end

            % Bottom neighbor (next row, same column) - add to lower triangle
            if i < m
                neighbor = i * m + j;
                ent = ent + 1;
                rows(ent) = neighbor;
                cols(ent) = node;
                vals(ent) = -1;
            end
        end
    end

    % Create sparse matrix from triplet format (lower triangle)
    A = sparse(rows(1:ent), cols(1:ent), vals(1:ent), n, n);

    % Symmetrize: A = L + L' where L is strictly lower triangular + diagonal
    % This creates the full matrix with correct off-diagonal values
    A = A + A' - spdiags(diag(A), 0, n, n);
end
