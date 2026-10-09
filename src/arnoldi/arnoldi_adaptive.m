function [V, lastv, H, h, beta, succ_flag, num_oracle, S] = ...
    arnoldi_adaptive(S, A, V, m, H, start_ind, param)
%ARNOLDI_ADAPTIVE Stop a truncated Arnoldi cycle using a sketch condition.
% The sketch monitors the basis; Arnoldi products use Euclidean inner products.
% succ_flag is false when a normalization breakdown is detected.

H(m + 1, m) = 0;
if size(V, 2) < m
    V(:, m) = 0;
end
trunc = param.truncation_length;
succ_flag = true;
num_oracle = 0;
n = size(V, 1);
s0 = 30;
if isempty(S)
    S = sketching_mat(s0, n, param.sketching_mat_type);
end
SV = zeros(size(S, 1), m + 1);

lastv = V(:, start_ind);
beta = norm(lastv);
if beta == 0
    succ_flag = false;
    h = 0;
    ell = start_ind - 1;
    H = H(1 : ell, 1 : ell);
    V = V(:, 1 : ell);
    return
end
lastv = lastv / beta;
V(:, start_ind) = lastv;
H(start_ind, 1 : start_ind - 1) = beta * H(start_ind, 1 : start_ind - 1);
SV(:, 1 : start_ind) = S * V(:, 1 : start_ind);
lastp = SV(:, start_ind);
examine_period = max(1, min(ceil((m - start_ind + 1) / 10), 5));

for j = start_ind : m
    if isnumeric(A)
        lastv = A * lastv;
    else
        lastv = A(lastv);
    end
    num_oracle = num_oracle + 1;

    j_trunc_start = max(1, j - trunc + 1);
    for i = j_trunc_start : j
        ip = V(:, i)' * lastv;
        H(i, j) = H(i, j) + ip;
        lastv = lastv - V(:, i) * ip;
    end

    H(j + 1, j) = norm(lastv);
    if H(j + 1, j) == 0 || ...
            abs(H(j + 1, j)) < j * eps * norm(H(1 : j + 1, j))
        succ_flag = false;
        break
    end
    lastv = lastv / H(j + 1, j);
    lastp = S * lastv;
    if j < m
        V(:, j + 1) = lastv;
        SV(:, j + 1) = lastp;
        if size(S, 1) < param.sketching_size_control * (j + 1)
            S_incr = sketching_mat(s0, n, param.sketching_mat_type);
            S = [S; S_incr];
            SV = [SV; S_incr * V(:, 1 : j + 1), ...
                zeros(size(S_incr, 1), size(SV, 2) - j - 1)];
            lastp = SV(:, j + 1);
            if param.verbose >= 2
                fprintf('sketching size increased to %d\n', size(S, 1));
            end
        end
        if mod(num_oracle, examine_period) == 0
            if cond(SV(:, 1 : j + 1)) > param.cond_tol
                break
            end
        end
    end
end

m_now = j;
h = H(m_now + 1, m_now);
H = H(1 : m_now, 1 : m_now);
V = V(:, 1 : m_now);
SV = SV(:, 1 : m_now);
if ~succ_flag
    return
end

switch param.update
    case 'none'
    case 'last_orth'
        [lastv, H, h] = arnoldi_last_orth_update(V, m_now, lastv, H, h);
    case 'last_sorth'
        [lastv, H, h, ~] = arnoldi_last_sorth_update(V, m_now, lastv, H, h, SV, lastp);
    case 'last_hmorth'
        [lastv, H, h] = arnoldi_last_hmorth_update(V, m_now, lastv, H, h);
    case 'last_shmorth'
        [lastv, H, h, ~] = arnoldi_last_shmorth_update(V, m_now, lastv, H, h, SV, lastp);
    otherwise
        error('FUNM_QUAD:arnoldi_adaptive', 'Unknown adaptive Arnoldi update.');
end
if h == 0 || abs(h) < m_now * eps * norm([H(:, m_now); h])
    succ_flag = false;
end

end
