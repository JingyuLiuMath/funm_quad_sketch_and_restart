exp_conv_diff_log_settings;

file_prefix = "conv_diff_n_" + string(n) + "_nu_" + string(nu);
matrix_file = matrix_data_path + file_prefix + "_matrix.mat";
reference_file = matrix_data_path + file_prefix + "_sol.mat";

A = two_dimensional_conv_diff(n, nu);

b0 = ones(N, 1) / sqrt(N);
b = A * b0;

fprintf('n: %d, N: %d, nu: %.1e\n', n, N, nu);
t_ref_start = tic;
[F, exitflag] = logm(full(speye(N) + A));
f_ex = F * b0;
t_ref = toc(t_ref_start);
assert(exitflag == 0, 'logm did not complete successfully.');
fprintf('Reference computation time: %.1e seconds\n', t_ref);

save(matrix_file, ...
    "A", "b", "n", "N", "nu");
save(reference_file, ...
    "f_ex", "t_ref", "exitflag");
