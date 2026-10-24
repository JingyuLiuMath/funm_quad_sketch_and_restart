exp_graph_laplacian_invsqrt_settings;

fprintf("Graph: %s\n", graph_name);

matrix_file = matrix_data_path + "frac_laplacian_matrix_" + graph_name + ".mat";
reference_file = matrix_data_path + "frac_laplacian_sol_" + graph_name + ".mat";

load(matrix_file, 'A');
N = size(A, 1);
rng(0, 'twister');
b0 = complex(rand(N, 1), rand(N, 1));
b = A * b0;
save(matrix_file, 'b', 'b0', '-append');

t_ref_start = tic;

% For the Gnutella data, zero columns span the kernel.
% Extract the nonsingular block without reordering A.
I = find(full(any(A, 1)));
B = A(I, I);
[F, sqrtm_residual] = sqrtm(full(B));
f_I = F * b0(I);
clear F;
f_ex = A(:, I) * (B \ f_I);

t_ref = toc(t_ref_start);
save(reference_file, 'f_ex', 'graph_name', 't_ref', 'sqrtm_residual');
fprintf('sqrtm residual on the nonzero-spectrum block: %.3e\n', sqrtm_residual);
fprintf('Reference saved to %s (%.3g seconds)\n', reference_file, t_ref);
