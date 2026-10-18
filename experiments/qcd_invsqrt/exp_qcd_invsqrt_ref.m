exp_qcd_invsqrt_settings;

Q = lattice_qcd(qcd_matrix_file);
source = load(qcd_matrix_file, 'source_matrix_file');
source_matrix_file = source.source_matrix_file;
N = size(Q, 1);

% A(v) = Q*(Q*v) and b = Q*e1/||Q*e1|| give the normalized sign action.
b0 = zeros(N, 1);
b0(1) = 1;
b = Q * b0;
normb = norm(b);
b = b / normb;

% Reuse the supplied sign(Q)*e1 reference with the same normalization.
reference = load(qcd_reference_file, ...
    'qcd_nonhermitian_exact', 'source_reference_file');
qcd_nonhermitian_exact = reference.qcd_nonhermitian_exact;
source_reference_file = reference.source_reference_file;
f_ex = qcd_nonhermitian_exact / normb;

fprintf('QCD matrix: N=%d, nnz(Q)=%d\n', N, nnz(Q));
fprintf('Right-hand side normalization: %.16e\n', normb);
fprintf('Supplied reference: %s\n', source_reference_file);

save(data_path + "qcd_matrix.mat", ...
    "Q", "b", "normb", "source_matrix_file", '-v7.3');
save(data_path + "qcd_sol.mat", ...
    "f_ex", "qcd_nonhermitian_exact", "normb", "source_reference_file");
