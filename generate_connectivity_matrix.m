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
    %
    % Description:
    %   Creates a sparse matrix representing a grid connectivity problem.
    %   - Diagonal elements: 4 (each node has 4 neighbors)
    %   - Off-diagonal elements: -1 (for connected neighbors)
    %   This is a classic 5-point Laplacian stencil matrix from FEM/FDM
    %
    % Reference: Exercise 2 from "Algebra matricial computacional"

    n = m * m;

    % Initialize sparse matrix
    A = sparse(n, n);

    % Build connectivity matrix
    for i = 1:m
        for j = 1:m
            % Linear index for node (i,j)
            node = (i - 1) * m + j;

            % Diagonal element: 4
            A(node, node) = 4;

            % Connect to right neighbor (j+1)
            if j < m
                neighbor = (i - 1) * m + (j + 1);
                A(node, neighbor) = -1;
                A(neighbor, node) = -1;
            end

            % Connect to bottom neighbor (i+1)
            if i < m
                neighbor = i * m + j;
                A(node, neighbor) = -1;
                A(neighbor, node) = -1;
            end
        end
    end

    % Ensure symmetric matrix (MATLAB sparse might double some entries)
    A = sparse(A + A' - diag(diag(A)));
end
