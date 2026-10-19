exp_conv_diff_log_settings;

for krylov_dim = krylov_dim_list
    mode_list = ["fixed", "adaptive"];
    for sfactor = sfactor_list
        for mode = mode_list
            for thick = thick_list
                result_list = struct([]);
                it_result = 0;
                if mode == "fixed"
                    for it_method = 1:numel(fixed_method_list)
                        method = fixed_method_list(it_method);
                        if thick && (method == "hmarnoldi" || method == "shmarnoldi")
                            continue
                        end
                        factor = sfactor;
                        if method == "arnoldi" || method == "hmarnoldi"
                            factor = 0;
                        end
                        file_name = "conv_diff_log_" + mode + "_" + method ...
                            + "_sfactor_" + string(factor) ...
                            + "_thick_" + string(double(thick)) ...
                            + "_krylov_dim_" + string(krylov_dim) + ".mat";
                        load(data_path + file_name, "result");
                        it_repeat = find(result.t_list == median(result.t_list), 1);
                        info = result.info(it_repeat);
                        it_result = it_result + 1;
                        result_list(it_result).method_name = method;
                        result_list(it_result).res_hist = info.err_history;
                        result_list(it_result).method_index = it_method;
                        result_list(it_result).truncation_index = 1;
                    end
                else
                    for it_method = 1:numel(adaptive_update_list)
                        method = adaptive_update_list(it_method);
                        factor = sfactor;
                        for it_trunc = 1:numel(truncation_length_list)
                            truncation_length = truncation_length_list(it_trunc);
                            file_name = "conv_diff_log_" + mode + "_" + method ...
                                + "_sfactor_" + string(factor) ...
                                + "_trunc_" + string(truncation_length) ...
                                + "_thick_" + string(double(thick)) ...
                                + "_krylov_dim_" + string(krylov_dim) + ".mat";
                            load(data_path + file_name, "result");
                            it_repeat = find(result.t_list == median(result.t_list), 1);
                            info = result.info(it_repeat);
                            it_result = it_result + 1;
                            result_list(it_result).method_name = method + ", trunc=" + string(truncation_length);
                            result_list(it_result).res_hist = info.err_history;
                            result_list(it_result).method_index = it_method;
                            result_list(it_result).truncation_index = it_trunc;
                        end
                    end
                end
                file_name = "conv_diff_log_res_" + mode + "_sfactor_" + string(sfactor) ...
                    + "_thick_" + string(double(thick)) ...
                    + "_krylov_dim_" + string(krylov_dim);
                figure_title = mode + ", sfactor=" + string(sfactor) + ", thick=" + string(double(thick)) ...
                    + ", krylov_dim=" + string(krylov_dim);
                plot_res(result_list, figure_path + file_name, figure_title);
            end
        end
    end
end
