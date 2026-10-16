exp_conv_diff_log_settings;

A = two_dimensional_conv_diff(n, nu);

b0 = ones(N, 1) / sqrt(N);
b = A * b0;

% funm_quad_log(A,b) evaluates g(A)b = log(I+A)*b0 at this RHS.
fprintf('n: %d, N: %d, nu: %.1e\n', n, N, nu);
t_ref_start = tic;
[F, exitflag] = logm(full(speye(N) + A));
f_ex = F * b0;
t_ref = toc(t_ref_start);
assert(exitflag == 0, 'logm did not complete successfully.');
fprintf('Reference computation time: %.1e seconds\n', t_ref);

save(data_path + "convdiff_matrix.mat", ...
    "A", "b", "nu");
save(data_path + "convdiff_sol.mat", ...
    "f_ex", "t_ref", "exitflag");
