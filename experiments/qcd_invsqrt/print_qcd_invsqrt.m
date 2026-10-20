exp_qcd_invsqrt_settings;

for krylov_dim = krylov_dim_list
    mode_list = ["fixed", "adaptive"];
    for sfactor = sfactor_list
        for truncation_length = truncation_length_list
            fprintf('\nkrylov_dim = %d, sfactor = %d, trunc = %d (fixed results reused; unsketched sfactor=0)\n', ...
                krylov_dim, sfactor, truncation_length);
            fprintf('%-9s %-14s %7s %7s %5s %12s %-36s %16s %14s %22s\n', ...
                'mode', 'method', 'sfactor', 'trunc', 'thick', 'time', 'reason', ...
                'iteration number', 'relative error', 'matrix-vector products');
            for mode = mode_list
                if mode == "fixed"
                    method_list = fixed_method_list;
                else
                    method_list = adaptive_update_list;
                end
                for thick = thick_list
                    for method = method_list
                        if mode == "fixed" && thick ...
                                && (method == "hmarnoldi" || method == "shmarnoldi")
                            continue
                        end
                        factor = sfactor;
                        if mode == "fixed" && (method == "arnoldi" || method == "hmarnoldi")
                            factor = 0;
                        end
                        file_name = "qcd_invsqrt_" + mode + "_" + method ...
                            + "_krylov_dim_" + string(krylov_dim) ...
                            + "_sfactor_" + string(factor) ...
                            + "_thick_" + string(double(thick));
                        if mode == "adaptive"
                            file_name = file_name + "_trunc_" + string(truncation_length);
                        end
                        file_name = file_name + ".mat";
                        load(data_path + file_name, "result");
                        it_repeat = find(result.t_list == median(result.t_list), 1);
                        info = result.info(it_repeat);
                        fprintf('%-9s %-14s %7d %7g %5d %12.5e %-36s %16d %14.5e %22d\n', ...
                            mode, method, result.sfactor, result.param.truncation_length, thick, ...
                            info.t_total, info.stop_condition, info.iter, ...
                            info.final_err, sum(info.num_oracle));
                    end
                end
            end
        end
    end
end
