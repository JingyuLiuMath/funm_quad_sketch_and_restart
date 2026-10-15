function [arnoldi_data, succ_flag, num_oracle] = ...
    shmarnoldi_fix(op_S, A, arnoldi_data, m, start_ind, param)

if ~isa(A, 'function_handle')
    op_A = @(x) A * x;
else
    op_A = @(x) A(x);
end

arnoldi_data.H(m + 1, m) = 0;
if size(arnoldi_data.V, 2) < m
    arnoldi_data.V(:, m) = 0;
end
if size(arnoldi_data.SV, 2) < m
    arnoldi_data.SV(:, m) = 0;
end
if size(arnoldi_data.Z, 2) < m
    arnoldi_data.Z(:, m) = 0;
end
if size(arnoldi_data.SZ, 2) < m
    arnoldi_data.SZ(:, m) = 0;
end
trunc = param.truncation_length;
succ_flag = true;
num_oracle = 0;
m_init = start_ind - 1;

arnoldi_data.beta = norm(arnoldi_data.sv);
if arnoldi_data.beta == 0
    succ_flag = false;
    arnoldi_data.h = 0;
    arnoldi_data.H = arnoldi_data.H(1 : m_init, 1 : m_init);
    arnoldi_data.V = arnoldi_data.V(:, 1 : m_init);
    arnoldi_data.SV = arnoldi_data.SV(:, 1 : m_init);
    arnoldi_data.Z = arnoldi_data.Z(:, 1 : m_init);
    arnoldi_data.SZ = arnoldi_data.SZ(:, 1 : m_init);
    return
end
arnoldi_data.sv = arnoldi_data.sv / arnoldi_data.beta;
arnoldi_data.v = arnoldi_data.v / arnoldi_data.beta;
arnoldi_data.SV(:, start_ind) = arnoldi_data.sv;
arnoldi_data.V(:, start_ind) = arnoldi_data.v;
arnoldi_data.H(start_ind, 1 : m_init) = arnoldi_data.beta * arnoldi_data.H(start_ind, 1 : m_init);

for j = start_ind : m
    arnoldi_data.v = op_A(arnoldi_data.v);
    arnoldi_data.sv = op_S(arnoldi_data.v);
    num_oracle = num_oracle + 1;

    rho = arnoldi_data.sv' * arnoldi_data.SV(:, j);
    if abs(rho) <= eps * norm(arnoldi_data.sv) * norm(arnoldi_data.SV(:, j))
        succ_flag = false;
        m = j;
        arnoldi_data.H(j + 1, j) = 0;
        break
    end
    arnoldi_data.Z(:, j) = arnoldi_data.v / conj(rho);
    arnoldi_data.SZ(:, j) = arnoldi_data.sv / conj(rho);

    j_trunc_start = max(1, j - trunc + 1);
    for i = j_trunc_start : j
        h_i_j = arnoldi_data.SZ(:, i)' * arnoldi_data.sv;
        arnoldi_data.H(i, j) = arnoldi_data.H(i, j) + h_i_j;
        arnoldi_data.sv = arnoldi_data.sv - arnoldi_data.SV(:, i) * h_i_j;
        arnoldi_data.v = arnoldi_data.v - arnoldi_data.V(:, i) * h_i_j;
    end

    arnoldi_data.H(j + 1, j) = norm(arnoldi_data.sv);
    if abs(arnoldi_data.H(j + 1, j)) <= j * eps * norm(arnoldi_data.H(1 : j + 1, j))
        succ_flag = false;
        m = j;
        break
    end
    arnoldi_data.sv = arnoldi_data.sv / arnoldi_data.H(j + 1, j);
    arnoldi_data.v = arnoldi_data.v / arnoldi_data.H(j + 1, j);
    if j < m
        arnoldi_data.SV(:, j + 1) = arnoldi_data.sv;
        arnoldi_data.V(:, j + 1) = arnoldi_data.v;
    end
end

arnoldi_data.h = arnoldi_data.H(m + 1, m);
arnoldi_data.H = arnoldi_data.H(1 : m, 1 : m);
arnoldi_data.V = arnoldi_data.V(:, 1 : m);
arnoldi_data.SV = arnoldi_data.SV(:, 1 : m);
arnoldi_data.Z = arnoldi_data.Z(:, 1 : m);
arnoldi_data.SZ = arnoldi_data.SZ(:, 1 : m);

end
