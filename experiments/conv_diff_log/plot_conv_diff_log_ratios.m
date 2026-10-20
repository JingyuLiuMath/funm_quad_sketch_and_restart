exp_conv_diff_log_settings;

for krylov_dim = krylov_dim_list
    load(data_path + "conv_diff_log_fixed_arnoldi_krylov_dim_" ...
        + string(krylov_dim) + "_sfactor_0_thick_0.mat", "result");
    it_repeat = find(result.t_list == median(result.t_list), 1);
    baseline_time = result.info(it_repeat).t_total;
    baseline_arnoldi_time = result.info(it_repeat).arnoldi_info.t_total ...
        / result.info(it_repeat).iter;

    mode_list = ["fixed", "adaptive"];
    for sfactor = sfactor_list
        for mode = mode_list
            if mode == "fixed"
                method_list = fixed_method_list;
                trunc_list = inf;
            else
                method_list = adaptive_update_list;
                trunc_list = truncation_length_list;
            end
            for thick = thick_list
                result_list = struct([]);
                it_method = 0;
                for truncation_length = trunc_list
                    for method = method_list
                        if mode == "fixed" && thick ...
                                && (method == "hmarnoldi" || method == "shmarnoldi")
                            continue
                        end
                        it_method = it_method + 1;
                        factor = sfactor;
                        if mode == "fixed" && (method == "arnoldi" || method == "hmarnoldi")
                            factor = 0;
                        end
                        file_name = "conv_diff_log_" + mode + "_" + method ...
                            + "_krylov_dim_" + string(krylov_dim) ...
                            + "_sfactor_" + string(factor) ...
                            + "_thick_" + string(double(thick));
                        if mode == "adaptive"
                            file_name = file_name + "_trunc_" + string(truncation_length);
                        end
                        file_name = file_name + ".mat";
                        load(data_path + file_name, "result");
                        it_repeat = find(result.t_list == median(result.t_list), 1);
                        result_list(it_method).method_name = method;
                        result_list(it_method).info = result.info(it_repeat);
                        if mode == "adaptive"
                            result_list(it_method).method_name = method + ", trunc=" + string(truncation_length);
                        end
                    end
                end
                file_name = "conv_diff_log_ratios_" + mode ...
                    + "_krylov_dim_" + string(krylov_dim) ...
                    + "_sfactor_" + string(sfactor) ...
                    + "_thick_" + string(double(thick));
                figure_prefix = figure_path + file_name;
                figure_title = "Convection-diffusion log, " + mode + ", sfactor=" + string(sfactor) ...
                    + ", thick=" + string(double(thick)) ...
                    + ", krylov_dim=" + string(krylov_dim);
                plot_ratios_total(result_list, figure_prefix, figure_title, baseline_time);
                if mode == "fixed"
                    plot_ratios_arnoldi(result_list, figure_prefix + "_arnoldi", ...
                        figure_title + ", average Arnoldi call", baseline_arnoldi_time);
                end
            end
        end
    end
end
