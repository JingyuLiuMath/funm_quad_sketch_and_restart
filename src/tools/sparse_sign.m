function S = sparse_sign(s, N)

zeta = min(s, 8);

nnzTotal = N * zeta;
rows = zeros(nnzTotal, 1);
cols_idx = zeros(nnzTotal, 1);
vals = zeros(nnzTotal, 1);

ptr = 1;
for j = 1:N
    ind_j = randperm(s, zeta);
    signs = sign(randn(zeta,1));
    global_ind = ptr:(ptr+zeta-1);
    rows(global_ind) = ind_j(:);
    cols_idx(global_ind) = j;
    vals(global_ind) = signs;
    ptr = ptr + zeta;
end

vals = vals / sqrt(zeta); 
S = sparse(rows, cols_idx, vals, s, N);

end