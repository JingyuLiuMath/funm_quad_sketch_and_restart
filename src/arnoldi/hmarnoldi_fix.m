function [V, lastv, H, h, beta, succ_flag, num_oracle] = ...
    hmarnoldi_fix(A, V, m, H, start_ind, param)
%HMARNOLDI_FIX Harmonic Arnoldi with full-space Euclidean pairings.

H(m + 1, m) = 0;
if size(V, 2) < m
    V(:, m) = 0;
end
trunc = param.truncation_length;
succ_flag = true;
num_oracle = 0;
ell = start_ind - 1;

Z = zeros(size(V, 1), m);
lastv = V(:, start_ind);
beta = norm(lastv);
if beta == 0
    succ_flag = false;
    h = 0;
    H = H(1 : ell, 1 : ell);
    V = V(:, 1 : ell);
    return
end
lastv = lastv / beta;
V(:, start_ind) = lastv;
H(start_ind, 1 : start_ind - 1) = beta * H(start_ind, 1 : start_ind - 1);
if ell > 0
    % Reconstruct A*Vprefix from the seeded relation, without A calls.
    % The joint dual satisfies Zprefix' * Vprefix = I.
    Q = V(:, 1 : start_ind) * H(1 : start_ind, 1 : ell);
    G = V(:, 1 : ell)' * Q;
    Z(:, 1 : ell) = Q / G;
end

for j = start_ind : m
    if isnumeric(A)
        lastv = A * lastv;
    else
        lastv = A(lastv);
    end
    num_oracle = num_oracle + 1;

    rho = lastv' * V(:, j);
    if rho == 0 || abs(rho) < eps * norm(lastv) * norm(V(:, j))
        succ_flag = false;
        m = j;
        H(j + 1, j) = 0;
        break
    end
    Z(:, j) = lastv / conj(rho);

    j_trunc_start = max(1, j - trunc + 1);
    for i = j_trunc_start : j
        ip = Z(:, i)' * lastv;
        H(i, j) = H(i, j) + ip;
        lastv = lastv - V(:, i) * ip;
    end

    H(j + 1, j) = norm(lastv);
    if H(j + 1, j) == 0 || abs(H(j + 1, j)) < j * eps * norm(H(1 : j + 1, j))
        succ_flag = false;
        m = j;
        break
    end
    lastv = lastv / H(j + 1, j);
    if j < m
        V(:, j + 1) = lastv;
    end
end

h = H(m + 1, m);
H = H(1 : m, 1 : m);
V = V(:, 1 : m);

end
