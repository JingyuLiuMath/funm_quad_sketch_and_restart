exp_graph_laplacian_invsqrt_settings;

fprintf("Graph: %s\n", graph_name);

matrix_file = matrix_data_path + "frac_laplacian_matrix_" + graph_name + ".mat";
reference_file = matrix_data_path + "frac_laplacian_sol_" + graph_name + ".mat";

load(matrix_file, 'A', 'matrix_name');
N = size(A, 1);
rng(0, 'twister');
b0 = complex(rand(N, 1), rand(N, 1));
b = A * b0;
save(matrix_file, 'A', 'b', 'b0', 'matrix_name', '-v7.3');
t_ref_start = tic;
[F, sqrtm_residual] = sqrtm(full(A));
f_ex = F * b0;
t_ref = toc(t_ref_start);
save(reference_file, 'f_ex', 'graph_name', 't_ref', 'sqrtm_residual');
fprintf('sqrtm relative matrix residual: %.3e\n', sqrtm_residual);
fprintf('Reference saved to %s (%.3g seconds)\n', reference_file, t_ref);
