exp_conv_diff_log_settings;

h = 1 / (n + 1);
N = n^2;

% Preserve the convection-diffusion operator and its ordering.
e = ones(n, 1);
L = spdiags([-e, 2 * e, -e], -1 : 1, n, n);
C = spdiags([-e, e], [-1, 0], n, n);
G = D / h^2 * L + C / h;
A = kron(speye(n), G) + kron(G.', speye(n));

b0 = ones(N, 1) / sqrt(N);
b = A * b0;

% funm_quad_log(A,b) evaluates g(A)b = log(I+A)*b0 at this RHS.
fprintf('n: %d, N: %d, D: %.1e\n', n, N, D);
t_ref_start = tic;
[F, exitflag] = logm(full(speye(N) + A));
f_ex = F * b0;
t_ref = toc(t_ref_start);
clear F
assert(exitflag == 0, 'logm did not complete successfully.');
fprintf('Reference computation time: %.2f seconds\n', t_ref);

save(data_path + "convdiff_matrix.mat", ...
    "A", "b", "D", "h");
save(data_path + "convdiff_sol.mat", ...
    "f_ex", "t_ref", "exitflag");
