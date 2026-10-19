exp_conv_diff_log_settings;

load(data_path + "convdiff_matrix.mat", "A", "b");
load(data_path + "convdiff_sol.mat", "f_ex");
param.exact = f_ex;
N = size(A, 1);

for krylov_dim = krylov_dim_list
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
                        [f, info] = funm_quad_log(A, b, myparam);
                        result.t_list(it_repeat) = info.t_total;
                        result.err_list(it_repeat) = info.final_err;
                        result.f(:, it_repeat) = f;
                        result.info(it_repeat) = info;

                        fprintf("adaptive %s, sfactor=%d, trunc=%d, thick=%d, krylov_dim=%d, run %d/%d\n", ...
                            method, sfactor, truncation_length, thick, krylov_dim, it_repeat, num_repeat);
                        fprintf("  time: %.1e\n", result.t_list(it_repeat));
                        fprintf("  reason: %s\n", info.stop_condition);
                        fprintf("  iteration number: %d\n", info.iter);
                        fprintf("  true relative error: %.1e\n", result.err_list(it_repeat));
                        fprintf("  matrix-vector products: %d\n", sum(info.num_oracle));
                    end

                    save(data_path + "conv_diff_log_adaptive_" + method ...
                        + "_sfactor_" + string(sfactor) ...
                        + "_trunc_" + string(truncation_length) ...
                        + "_thick_" + string(double(thick)) ...
                        + "_krylov_dim_" + string(krylov_dim) + ".mat", ...
                        "result", '-v7.3');
                end
            end
        end
    end
end
