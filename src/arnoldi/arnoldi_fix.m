function [arnoldi_data, succ_flag, num_oracle, arnoldi_info] = ...
    arnoldi_fix(A, arnoldi_data, m, start_ind, param)

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
arnoldi_data.H(m + 1, m) = 0;
if size(arnoldi_data.V, 2) < m
    arnoldi_data.V(:, m) = 0;
end
trunc = param.truncation_length;
succ_flag = true;
num_oracle = 0;
m_init = start_ind - 1;

arnoldi_data.beta = norm(arnoldi_data.v);
if arnoldi_data.beta == 0
    succ_flag = false;
    arnoldi_data.h = 0;
    arnoldi_data.H = arnoldi_data.H(1 : m_init, 1 : m_init);
    arnoldi_data.V = arnoldi_data.V(:, 1 : m_init);
    arnoldi_info.t_total = toc(t_total_start);
    arnoldi_info.t_other = arnoldi_info.t_total - ...
        (arnoldi_info.t_apply_op + arnoldi_info.t_apply_sketch + ...
        arnoldi_info.t_coeff + arnoldi_info.t_substract + arnoldi_info.t_update);
    return
end
arnoldi_data.v = arnoldi_data.v / arnoldi_data.beta;
arnoldi_data.V(:, start_ind) = arnoldi_data.v;
arnoldi_data.H(start_ind, 1 : m_init) = arnoldi_data.beta * arnoldi_data.H(start_ind, 1 : m_init);

for j = start_ind : m
    t_apply_op_start = tic;
    arnoldi_data.v = op_A(arnoldi_data.v);
    arnoldi_info.t_apply_op = arnoldi_info.t_apply_op + toc(t_apply_op_start);
    num_oracle = num_oracle + 1;

    for i = max(1, j - trunc + 1) : j
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
        m = j;
        break
    end
    
    arnoldi_data.v = arnoldi_data.v / arnoldi_data.H(j + 1, j);
    if j < m
        arnoldi_data.V(:, j + 1) = arnoldi_data.v;
    end
end

arnoldi_data.h = arnoldi_data.H(m + 1, m);
arnoldi_data.H = arnoldi_data.H(1 : m, 1 : m);
arnoldi_data.V = arnoldi_data.V(:, 1 : m);

arnoldi_info.t_total = toc(t_total_start);
arnoldi_info.t_other = arnoldi_info.t_total - ...
    (arnoldi_info.t_apply_op + arnoldi_info.t_apply_sketch + ...
    arnoldi_info.t_coeff + arnoldi_info.t_substract + arnoldi_info.t_update);

end
