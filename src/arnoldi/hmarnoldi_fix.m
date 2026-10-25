function [beta, v, V, h, H, succ_flag, num_oracle, arnoldi_info] = ...
    hmarnoldi_fix(A, v, m, param)

% Harmonic Arnoldi always starts from the first basis vector.

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
H = zeros(m, m, 'like', v);
Z = zeros(size(v, 1), m, 'like', v);
trunc = param.truncation_length;
succ_flag = true;
num_oracle = 0;

beta = norm(v);
if beta == 0
    succ_flag = false;
    h = 0;
    H = H([], []);
    V = V(:, []);
    arnoldi_info.t_total = toc(t_total_start);
    arnoldi_info.t_other = arnoldi_info.t_total - ...
        (arnoldi_info.t_apply_op + arnoldi_info.t_apply_sketch + ...
        arnoldi_info.t_coeff + arnoldi_info.t_substract + arnoldi_info.t_update);
    return
end
v = v / beta;
V(:, 1) = v;

for j = 1 : m
    t_apply_op_start = tic;
    v = op_A(v);
    arnoldi_info.t_apply_op = arnoldi_info.t_apply_op + toc(t_apply_op_start);
    num_oracle = num_oracle + 1;

    rho = v' * V(:, j);
    if abs(rho) <= eps * norm(v) * norm(V(:, j))
        succ_flag = false;
        h = 0;
        m = j;
        break
    end
    Z(:, j) = v / conj(rho);

    for i = max(1, j - trunc + 1) : j
        t_coeff_start = tic;
        h_i_j = Z(:, i)' * v;
        H(i, j) = H(i, j) + h_i_j;
        arnoldi_info.t_coeff = arnoldi_info.t_coeff + toc(t_coeff_start);
        t_substract_start = tic;
        v = v - V(:, i) * h_i_j;
        arnoldi_info.t_substract = arnoldi_info.t_substract + toc(t_substract_start);
    end

    h = norm(v);
    if abs(h) <= j * eps * norm([H(1 : j, j); h])
        succ_flag = false;
        m = j;
        break
    end

    v = v / h;
    if j < m
        H(j + 1, j) = h;
        V(:, j + 1) = v;
    end
end

if ~succ_flag && m < size(H, 1)
    H = H(1 : m, 1 : m);
    V = V(:, 1 : m);
end

arnoldi_info.t_total = toc(t_total_start);
arnoldi_info.t_other = arnoldi_info.t_total - ...
    (arnoldi_info.t_apply_op + arnoldi_info.t_apply_sketch + ...
    arnoldi_info.t_coeff + arnoldi_info.t_substract + arnoldi_info.t_update);

end
