function [S_op, m_now, arnoldi_data, succ_flag, num_oracle, num_resketch, arnoldi_info] = ...
    arnoldi_adaptive(S_op, A, arnoldi_data, m, start_ind, param)

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
N = size(arnoldi_data.v, 1);
s = size(arnoldi_data.sv, 1);
arnoldi_data.H(m + 1, m) = 0;
if size(arnoldi_data.V, 2) < m
    arnoldi_data.V(:, m) = 0;
end
if size(arnoldi_data.SV, 2) < m
    arnoldi_data.SV(:, m) = 0;
end
trunc = param.truncation_length;
succ_flag = true;
num_oracle = 0;
num_resketch = 0;
m_init = start_ind - 1;

arnoldi_data.beta = norm(arnoldi_data.v);
if arnoldi_data.beta == 0
    succ_flag = false;
    m_now = m_init;
    arnoldi_data.h = 0;
    arnoldi_data.H = arnoldi_data.H(1 : m_init, 1 : m_init);
    arnoldi_data.V = arnoldi_data.V(:, 1 : m_init);
    arnoldi_data.SV = arnoldi_data.SV(:, 1 : m_init);
    arnoldi_info.t_total = toc(t_total_start);
    arnoldi_info.t_other = arnoldi_info.t_total;
    return
end
arnoldi_data.v = arnoldi_data.v / arnoldi_data.beta;
arnoldi_data.sv = arnoldi_data.sv / arnoldi_data.beta;
arnoldi_data.V(:, start_ind) = arnoldi_data.v;
arnoldi_data.SV(:, start_ind) = arnoldi_data.sv;
arnoldi_data.H(start_ind, 1 : m_init) = arnoldi_data.beta * arnoldi_data.H(start_ind, 1 : m_init);
examine_period = 1;

for j = start_ind : m
    t_apply_op_start = tic;
    arnoldi_data.v = op_A(arnoldi_data.v);
    arnoldi_info.t_apply_op = arnoldi_info.t_apply_op + toc(t_apply_op_start);
    num_oracle = num_oracle + 1;

    j_trunc_start = max(1, j - trunc + 1);
    for i = j_trunc_start : j
        t_coeff_start = tic;
        h_i_j = arnoldi_data.V(:, i)' * arnoldi_data.v;
        arnoldi_data.H(i, j) = arnoldi_data.H(i, j) + h_i_j;
        arnoldi_info.t_coeff = arnoldi_info.t_coeff + toc(t_coeff_start);
        t_substract_start = tic;
        arnoldi_data.v = arnoldi_data.v - arnoldi_data.V(:, i) * h_i_j;
        arnoldi_info.t_substract = arnoldi_info.t_substract + toc(t_substract_start);
    end

    arnoldi_data.H(j + 1, j) = norm(arnoldi_data.v);
    if abs(arnoldi_data.H(j + 1, j)) <= j * eps * norm(arnoldi_data.H(1 : j + 1, j))
        succ_flag = false;
        t_apply_sketch_start = tic;
        arnoldi_data.sv = S_op(arnoldi_data.v);
        arnoldi_info.t_apply_sketch = arnoldi_info.t_apply_sketch + toc(t_apply_sketch_start);
        break
    end

    arnoldi_data.v = arnoldi_data.v / arnoldi_data.H(j + 1, j);
    t_apply_sketch_start = tic;
    arnoldi_data.sv = S_op(arnoldi_data.v);
    arnoldi_info.t_apply_sketch = arnoldi_info.t_apply_sketch + toc(t_apply_sketch_start);
    if j < m
        arnoldi_data.V(:, j + 1) = arnoldi_data.v;
        arnoldi_data.SV(:, j + 1) = arnoldi_data.sv;
        if s < param.sketching_size_control * (j + 1)
            s = s + s0;
            S_op = sparse_sign(s, N);
            num_resketch = num_resketch + 1;
            arnoldi_data.SV = zeros(s, m);
            t_apply_sketch_start = tic;
            arnoldi_data.SV(:, 1 : (j + 1)) = S_op(arnoldi_data.V(:, 1 : (j + 1)));
            arnoldi_info.t_apply_sketch = arnoldi_info.t_apply_sketch + toc(t_apply_sketch_start);
            arnoldi_data.sv = arnoldi_data.SV(:, j + 1);
        end
        if mod(num_oracle, examine_period) == 0
            if cond(arnoldi_data.SV(:, 1 : j + 1)) > param.cond_tol
                break
            end
        end
    end
end

m_now = j;
arnoldi_data.h = arnoldi_data.H(m_now + 1, m_now);
arnoldi_data.H = arnoldi_data.H(1 : m_now, 1 : m_now);
arnoldi_data.V = arnoldi_data.V(:, 1 : m_now);
arnoldi_data.SV = arnoldi_data.SV(:, 1 : m_now);
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
        arnoldi_data = arnoldi_last_orth_update(arnoldi_data);
    case "last_sorth"
        update_flag = true;
        t_update_start = tic;
        arnoldi_data = arnoldi_last_sorth_update(arnoldi_data);
    case "last_hmorth"
        update_flag = true;
        t_update_start = tic;
        arnoldi_data = arnoldi_last_hmorth_update(arnoldi_data);
    case "last_shmorth"
        update_flag = true;
        t_update_start = tic;
        arnoldi_data = arnoldi_last_shmorth_update(arnoldi_data);
end

if update_flag
    if abs(arnoldi_data.h) <= m_now * eps * norm([arnoldi_data.H(:, m_now); arnoldi_data.h])
        succ_flag = false;
    end
    arnoldi_info.t_update = toc(t_update_start);
end

arnoldi_info.t_total = toc(t_total_start);
arnoldi_info.t_other = arnoldi_info.t_total - ...
    (arnoldi_info.t_apply_op + arnoldi_info.t_apply_sketch + ...
    arnoldi_info.t_coeff + arnoldi_info.t_substract + arnoldi_info.t_update);

end
