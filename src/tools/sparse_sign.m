function op_S = sparse_sign(s, N)

zeta = 8;

nnz_S = N * zeta;
rows = zeros(nnz_S, 1);
cols_idx = zeros(nnz_S, 1);
vals = zeros(nnz_S, 1);

ptr = 1;
for j = 1 : N
    ind_j = randperm(s, zeta);
    signs = sign(randn(zeta, 1));
    global_ind = ptr : (ptr + zeta - 1);
    rows(global_ind) = ind_j(:);
    cols_idx(global_ind) = j;
    vals(global_ind) = signs;
    ptr = ptr + zeta;
end

vals = vals / sqrt(zeta); 
S = sparse(rows, cols_idx, vals, s, N);
op_S = @(x) S * x;

end