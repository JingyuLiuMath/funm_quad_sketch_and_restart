function [beta, v, V, sv, SV, h, H, succ_flag, num_oracle, arnoldi_info] = ...
    sarnoldi_fix(op_S, A, v, V_init, sv, SV_init, H_init, m, start_ind, param)

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
V = zeros(size(v, 1), m, 'like', v);
if ~isempty(V_init)
    V(:, 1 : size(V_init, 2)) = V_init;
end
H = zeros(m, m, 'like', v);
if ~isempty(H_init)
    H(1 : size(H_init, 1), 1 : size(H_init, 2)) = H_init;
end
SV = zeros(size(sv, 1), m, 'like', sv);
if ~isempty(SV_init)
    SV(:, 1 : size(SV_init, 2)) = SV_init;
end
trunc = param.truncation_length;
succ_flag = true;
num_oracle = 0;
m_init = start_ind - 1;

beta = norm(sv);
if beta == 0
    succ_flag = false;
    h = 0;
    H = H(1 : m_init, 1 : m_init);
    V = V(:, 1 : m_init);
    SV = SV(:, 1 : m_init);
    arnoldi_info.t_total = toc(t_total_start);
    arnoldi_info.t_other = arnoldi_info.t_total - ...
        (arnoldi_info.t_apply_op + arnoldi_info.t_apply_sketch + ...
        arnoldi_info.t_coeff + arnoldi_info.t_substract + arnoldi_info.t_update);
    return
end
sv = sv / beta;
v = v / beta;
SV(:, start_ind) = sv;
V(:, start_ind) = v;
H(start_ind, 1 : m_init) = beta * H(start_ind, 1 : m_init);

for j = start_ind : m
    t_apply_op_start = tic;
    v = op_A(v);
    arnoldi_info.t_apply_op = arnoldi_info.t_apply_op + toc(t_apply_op_start);
    t_apply_sketch_start = tic;
    sv = op_S(v);
    arnoldi_info.t_apply_sketch = arnoldi_info.t_apply_sketch + toc(t_apply_sketch_start);
    num_oracle = num_oracle + 1;

    j_trunc_start = max(1, j - trunc + 1);
    for i = j_trunc_start : j
        t_coeff_start = tic;
        h_i_j = SV(:, i)' * sv;
        H(i, j) = H(i, j) + h_i_j;
        arnoldi_info.t_coeff = arnoldi_info.t_coeff + toc(t_coeff_start);
        t_substract_start = tic;
        sv = sv - SV(:, i) * h_i_j;
        v = v - V(:, i) * h_i_j;
        arnoldi_info.t_substract = arnoldi_info.t_substract + toc(t_substract_start);
    end

    h = norm(sv);
    if abs(h) <= j * eps * norm([H(1 : j, j); h])
        succ_flag = false;
        m = j;
        break
    end

    sv = sv / h;
    v = v / h;
    if j < m
        H(j + 1, j) = h;
        SV(:, j + 1) = sv;
        V(:, j + 1) = v;
    end
end

if ~succ_flag && m < size(H, 1)
    H = H(1 : m, 1 : m);
    V = V(:, 1 : m);
    SV = SV(:, 1 : m);
end

arnoldi_info.t_total = toc(t_total_start);
arnoldi_info.t_other = arnoldi_info.t_total - ...
    (arnoldi_info.t_apply_op + arnoldi_info.t_apply_sketch + ...
    arnoldi_info.t_coeff + arnoldi_info.t_substract + arnoldi_info.t_update);

end
