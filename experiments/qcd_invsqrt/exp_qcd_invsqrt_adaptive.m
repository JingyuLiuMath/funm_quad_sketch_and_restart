exp_qcd_invsqrt_settings;

matrix_file = matrix_data_path + "qcd_matrix.mat";
reference_file = matrix_data_path + "qcd_sol.mat";

load(matrix_file, "Q", "b");
load(reference_file, "f_ex");
param.exact = f_ex;
N = size(Q, 1);
A = @(v) Q * (Q * v);
norm_exact = norm(f_ex);

krylov_dim = adaptive_krylov_dim;
for sfactor = sfactor_list
    for truncation_length = truncation_length_list
        for method = adaptive_update_list
            for thick = thick_list
                myparam = param;
                myparam.mode = "adaptive";
                myparam.krylov_dim = krylov_dim;
                myparam.arnoldi = "arnoldi";
                myparam.update = method;
                myparam.truncation_length = truncation_length;
                myparam.sketching_size_control = sfactor;
                myparam.thick = thick;

                result = struct();
                result.N = N;
                result.method = method;
                result.sfactor = sfactor;
                result.krylov_dim = krylov_dim;
                result.param = myparam;
                result.seed = seed_list;
                result.t_list = zeros(1, num_repeat);
                result.err_list = zeros(1, num_repeat);
                result.f = zeros(N, num_repeat);

                for it_repeat = 1 : num_repeat
                    rng(seed_list(it_repeat), 'twister');
                    [f, info] = funm_quad_invsqrt(A, b, myparam);
                    result.t_list(it_repeat) = info.t_total;
                    result.err_list(it_repeat) = norm(f - f_ex) / norm_exact;
                    result.f(:, it_repeat) = f;
                    result.info(it_repeat) = info;

                    fprintf("adaptive %s, sfactor=%d, trunc=%d, thick=%d, krylov_dim=%d, run %d/%d\n", ...
                        method, sfactor, truncation_length, thick, krylov_dim, it_repeat, num_repeat);
                    fprintf('  time: %.1e\n', info.t_total);
                    fprintf('  reason: %s\n', info.stop_condition);
                    fprintf('  iteration number: %d\n', info.iter);
                    fprintf('  relative error: %.1e\n', info.final_err);
                    fprintf('  matrix-vector products: %d\n', sum(info.num_oracle));
                end

                save(data_path + "qcd_invsqrt_adaptive_" + method ...
                    + "_krylov_dim_" + string(krylov_dim) ...
                    + "_sfactor_" + string(sfactor) ...
                    + "_thick_" + string(double(thick)) ...
                    + "_trunc_" + string(truncation_length) + ".mat", ...
                    "result", '-v7.3');
            end
        end
    end
end
