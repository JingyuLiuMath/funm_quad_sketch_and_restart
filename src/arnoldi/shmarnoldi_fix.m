function [beta, v, V, sv, SV, h, H, succ_flag, num_oracle, arnoldi_info] = ...
    shmarnoldi_fix(op_S, A, v, sv, SV_init, m, param)

% Sketched harmonic Arnoldi always starts from the first basis vector.

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
SV = zeros(size(sv, 1), m, 'like', sv);
if ~isempty(SV_init)
    SV(:, 1 : size(SV_init, 2)) = SV_init;
end
% Z = zeros(size(V, 1), m, 'like', v);
SZ = zeros(size(SV, 1), m, 'like', sv);
trunc = param.truncation_length;
succ_flag = true;
num_oracle = 0;

beta = norm(sv);
if beta == 0
    succ_flag = false;
    h = 0;
    H = H([], []);
    V = V(:, []);
    SV = SV(:, []);
    % Z = Z(:, []);
    % SZ = SZ(:, []);
    arnoldi_info.t_total = toc(t_total_start);
    arnoldi_info.t_other = arnoldi_info.t_total - ...
        (arnoldi_info.t_apply_op + arnoldi_info.t_apply_sketch + ...
        arnoldi_info.t_coeff + arnoldi_info.t_substract + arnoldi_info.t_update);
    return
end
sv = sv / beta;
v = v / beta;
SV(:, 1) = sv;
V(:, 1) = v;

for j = 1 : m
    t_apply_op_start = tic;
    v = op_A(v);
    arnoldi_info.t_apply_op = arnoldi_info.t_apply_op + toc(t_apply_op_start);
    t_apply_sketch_start = tic;
    sv = op_S(v);
    arnoldi_info.t_apply_sketch = arnoldi_info.t_apply_sketch + toc(t_apply_sketch_start);
    num_oracle = num_oracle + 1;

    rho = sv' * SV(:, j);
    if abs(rho) <= eps * norm(sv) * norm(SV(:, j))
        succ_flag = false;
        h = 0;
        m = j;
        break
    end
    % Z(:, j) = v / conj(rho);
    SZ(:, j) = sv / conj(rho);

    j_trunc_start = max(1, j - trunc + 1);
    for i = j_trunc_start : j
        t_coeff_start = tic;
        h_i_j = SZ(:, i)' * sv;
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
