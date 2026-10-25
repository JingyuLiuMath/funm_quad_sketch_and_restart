function [beta, v, V, sv, SV, h, H, op_S, m_now, ...
    succ_flag, num_oracle, num_resketch, arnoldi_info] = ...
    arnoldi_adaptive(op_S, A, v, V_init, sv, SV_init, H_init, m, start_ind, param)

% V_init, H_init, and SV_init contain the retained restart data.

if ~isa(A, 'function_handle')
    op_A = @(x) A * x;
else
    op_A = @(x) A(x);
end

arnoldi_info.t_total = 0;
arnoldi_info.t_apply_op = 0;
arnoldi_info.t_apply_sketch = 0;
arnoldi_info.t_coeff = 0;
arnoldi_info.t_substract = 0;
arnoldi_info.t_update = 0;
arnoldi_info.t_other = 0;

t_total_start = tic;
s0 = 30;
N = size(v, 1);
s = size(sv, 1);
H = zeros(m + 1, m, 'like', v);
if ~isempty(H_init)
    H(1 : size(H_init, 1), 1 : size(H_init, 2)) = H_init;
end
H(m + 1, m) = 0;
V = zeros(N, m, 'like', v);
if ~isempty(V_init)
    V(:, 1 : size(V_init, 2)) = V_init;
end
SV = zeros(s, m, 'like', sv);
if ~isempty(SV_init)
    SV(:, 1 : size(SV_init, 2)) = SV_init;
end
trunc = param.truncation_length;
succ_flag = true;
num_oracle = 0;
num_resketch = 0;
m_init = start_ind - 1;

beta = norm(v);
if beta == 0
    succ_flag = false;
    m_now = m_init;
    h = 0;
    H = H(1 : m_init, 1 : m_init);
    V = V(:, 1 : m_init);
    SV = SV(:, 1 : m_init);
    arnoldi_info.t_total = toc(t_total_start);
    arnoldi_info.t_other = arnoldi_info.t_total;
    return
end
v = v / beta;
sv = sv / beta;
V(:, start_ind) = v;
SV(:, start_ind) = sv;
H(start_ind, 1 : m_init) = beta * H(start_ind, 1 : m_init);
examine_period = 1;

for j = start_ind : m
    t_apply_op_start = tic;
    v = op_A(v);
    arnoldi_info.t_apply_op = arnoldi_info.t_apply_op + toc(t_apply_op_start);
    num_oracle = num_oracle + 1;

    j_trunc_start = max(1, j - trunc + 1);
    for i = j_trunc_start : j
        t_coeff_start = tic;
        h_i_j = V(:, i)' * v;
        H(i, j) = H(i, j) + h_i_j;
        arnoldi_info.t_coeff = arnoldi_info.t_coeff + toc(t_coeff_start);
        t_substract_start = tic;
        v = v - V(:, i) * h_i_j;
        arnoldi_info.t_substract = arnoldi_info.t_substract + toc(t_substract_start);
    end

    H(j + 1, j) = norm(v);
    if abs(H(j + 1, j)) <= j * eps * norm(H(1 : j + 1, j))
        succ_flag = false;
        t_apply_sketch_start = tic;
        sv = op_S(v);
        arnoldi_info.t_apply_sketch = arnoldi_info.t_apply_sketch + toc(t_apply_sketch_start);
        break
    end

    v = v / H(j + 1, j);
    t_apply_sketch_start = tic;
    sv = op_S(v);
    arnoldi_info.t_apply_sketch = arnoldi_info.t_apply_sketch + toc(t_apply_sketch_start);
    if j < m
        V(:, j + 1) = v;
        SV(:, j + 1) = sv;
        if s < param.sketching_size_control * (j + 1)
            s = s + s0;
            op_S = sparse_sign(s, N);
            num_resketch = num_resketch + 1;
            SV = zeros(s, m, 'like', sv);
            t_apply_sketch_start = tic;
            SV(:, 1 : (j + 1)) = op_S(V(:, 1 : (j + 1)));
            arnoldi_info.t_apply_sketch = arnoldi_info.t_apply_sketch + toc(t_apply_sketch_start);
            sv = SV(:, j + 1);
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
if size(V, 2) > m_now
    V = V(:, 1 : m_now);
end
if size(SV, 2) > m_now
    SV = SV(:, 1 : m_now);
end
if ~succ_flag
    arnoldi_info.t_total = toc(t_total_start);
    arnoldi_info.t_other = arnoldi_info.t_total - ...
        (arnoldi_info.t_apply_op + arnoldi_info.t_apply_sketch + ...
        arnoldi_info.t_coeff + arnoldi_info.t_substract + arnoldi_info.t_update);
    return
end

update_flag = false;
switch param.update
    case "last_orth"
        update_flag = true;
        t_update_start = tic;
        [v, sv, h, H] = arnoldi_last_orth_update(v, h, V, H, sv, SV);
    case "last_sorth"
        update_flag = true;
        t_update_start = tic;
        [v, sv, h, H] = arnoldi_last_sorth_update(v, h, V, H, sv, SV);
    case "last_hmorth"
        update_flag = true;
        t_update_start = tic;
        [v, sv, h, H] = arnoldi_last_hmorth_update(v, h, V, H, sv, SV);
    case "last_shmorth"
        update_flag = true;
        t_update_start = tic;
        [v, sv, h, H] = arnoldi_last_shmorth_update(v, h, V, H, sv, SV);
end

if update_flag
    if abs(h) <= m_now * eps * norm([H(:, m_now); h])
        succ_flag = false;
    end
    arnoldi_info.t_update = toc(t_update_start);
end

arnoldi_info.t_total = toc(t_total_start);
arnoldi_info.t_other = arnoldi_info.t_total - ...
    (arnoldi_info.t_apply_op + arnoldi_info.t_apply_sketch + ...
    arnoldi_info.t_coeff + arnoldi_info.t_substract + arnoldi_info.t_update);

end
