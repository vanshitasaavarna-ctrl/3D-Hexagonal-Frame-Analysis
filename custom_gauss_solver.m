% =========================================================================
% CVL7411 ADVANCED STRUCTURAL ANALYSIS - PROGRAMMING ASSIGNMENT GROUP (D)
% 3D SINGLE-STOREY HEXAGONAL SPACE FRAME WITH RIGID SLAB
% =========================================================================

function x = custom_gauss_solver(A, b)
n = length(b);
Aug = [A, b];

% Forward Elimination with Partial Pivoting
for k = 1:n-1
    [max_val, pivot_idx] = max(abs(Aug(k:n, k)));
    pivot_row = pivot_idx + k - 1;

    if max_val < 1e-12
        error('System matrix is singular or ill-conditioned.');
    end

    if pivot_row ~= k
        Aug([k, pivot_row], :) = Aug([pivot_row, k], :);
    end

    for i = k+1:n
        factor = Aug(i, k) / Aug(k, k);
        Aug(i, k:n+1) = Aug(i, k:n+1) - factor * Aug(k, k:n+1);
    end
end

% Back Substitution
x = zeros(n, 1);
x(n) = Aug(n, n+1) / Aug(n, n);
for i = n-1:-1:1
    x(i) = (Aug(i, n+1) - Aug(i, i+1:n) * x(i+1:n)) / Aug(i, i);
end
end