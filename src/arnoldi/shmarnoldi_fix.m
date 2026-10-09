function [V, lastv, H, h, beta, succ_flag, num_oracle] = shmarnoldi_fix(S, A, V, m, H, start_ind, param)
%SHMARNOLDI_FIX Sketched harmonic Arnoldi with a local basis.

H(m + 1, m) = 0;
if size(V, 2) < m
    V(:, m) = 0;
end
trunc = param.truncation_length;
succ_flag = true;
num_oracle = 0;
ell = start_ind - 1;

SV = zeros(size(S, 1), m);
Z = zeros(size(S, 1), m);
lastv = V(:, start_ind);
lastp = S * lastv;
beta = norm(lastp);
if beta == 0
    succ_flag = false;
    h = 0;
    H = H(1 : ell, 1 : ell);
    V = V(:, 1 : ell);
    return
end
lastp = lastp / beta;
lastv = lastv / beta;
V(:, start_ind) = lastv;
H(start_ind, 1 : start_ind - 1) = beta * H(start_ind, 1 : start_ind - 1);
SV(:, start_ind) = lastp;
if ell > 0
    SV(:, 1 : ell) = S * V(:, 1 : ell);
    % Reconstruct S*A*Vprefix from the seeded relation, without A calls.
    % The joint dual satisfies Zprefix' * SVprefix = I.
    Q = SV(:, 1 : start_ind) * H(1 : start_ind, 1 : ell);
    G = SV(:, 1 : ell)' * Q;
    Z(:, 1 : ell) = Q / G;
end

for j = start_ind : m
    if isnumeric(A)
        lastv = A * lastv;
    else
        lastv = A(lastv);
    end
    lastp = S * lastv;
    num_oracle = num_oracle + 1;

    rho = lastp' * SV(:, j);
    if rho == 0 || abs(rho) < eps * norm(lastp) * norm(SV(:, j))
        succ_flag = false;
        m = j;
        H(j + 1, j) = 0;
        break
    end
    Z(:, j) = lastp / conj(rho);

    j_trunc_start = max(1, j - trunc + 1);
    for i = j_trunc_start : j
        ip = Z(:, i)' * lastp;
        H(i, j) = H(i, j) + ip;
        lastp = lastp - SV(:, i) * ip;
        lastv = lastv - V(:, i) * ip;
    end

    H(j + 1, j) = norm(lastp);
    if H(j + 1, j) == 0 || abs(H(j + 1, j)) < j * eps * norm(H(1 : j + 1, j))
        succ_flag = false;
        m = j;
        break
    end
    lastp = lastp / H(j + 1, j);
    lastv = lastv / H(j + 1, j);
    if j < m
        SV(:, j + 1) = lastp;
        V(:, j + 1) = lastv;
    end
end

h = H(m + 1, m);
H = H(1 : m, 1 : m);
V = V(:, 1 : m);

end
