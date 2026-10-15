exp_conv_diff_log_settings;

load(data_path + "convdiff_matrix.mat", "A", "b");
load(data_path + "convdiff_sol.mat", "f_ex");
N = size(A, 1);
norm_f_ex = norm(f_ex);
param.exact = f_ex;

for method = fixed_method_list
    if method == "arnoldi" || method == "hmarnoldi"
        factor_list = 0; % No sketch: sfactor=0 denotes the baseline.
    else
        factor_list = sfactor_list;
    end
    for sfactor = factor_list
        myparam = param;
        myparam.mode = "fixed";
        myparam.truncation_length = inf;
        myparam.arnoldi = method;
        myparam.sketch_dim = sfactor * param.krylov_dim;
        for thick = thick_list
            if thick && (method == "hmarnoldi" || method == "shmarnoldi")
                continue
            end
            myparam.thick = thick;

            result = struct();
            result.n = N;
            result.method = method;
            result.sfactor = sfactor;
            result.param = myparam;
            result.seed = seed_list;
            result.t_list = zeros(1, num_repeat);
            result.err_list = zeros(1, num_repeat);
            result.f = zeros(N, num_repeat);

            for it_repeat = 1 : num_repeat
                rng(seed_list(it_repeat), 'twister');
                [f, info] = funm_quad_log(A, b, myparam);
                result.t_list(it_repeat) = info.t_total;
                result.err_list(it_repeat) = norm(f - f_ex) / norm_f_ex;
                result.f(:, it_repeat) = f;
                result.info(it_repeat) = info;

                fprintf("%s, sfactor=%d, thick=%d, run %d/%d\n", ...
                    method, sfactor, thick, it_repeat, num_repeat);
                fprintf("  time: %.1e\n", info.t_total);
                fprintf("  reason: %s\n", info.stop_condition);
                fprintf("  iteration number: %d\n", info.iter);
                fprintf("  true relative error: %.1e\n", result.err_list(it_repeat));
                fprintf("  matrix-vector products: %d\n", sum(info.num_oracle));
            end

            save(data_path + "conv_diff_log_fixed_" + method ...
                + "_sfactor_" + string(sfactor) ...
                + "_thick_" + string(double(thick)) + ".mat", ...
                "result", '-v7.3');
        end
    end
end
