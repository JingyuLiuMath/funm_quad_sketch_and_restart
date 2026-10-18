exp_qcd_invsqrt_settings;

load(data_path + "qcd_invsqrt_fixed_arnoldi_sfactor_0_thick_0.mat", "result");
it_repeat = find(result.t_list == median(result.t_list), 1);
baseline_time = result.info(it_repeat).t_total;

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
                    file_name = "qcd_invsqrt_" + mode + "_" + method ...
                        + "_sfactor_" + string(factor);
                    if mode == "adaptive"
                        file_name = file_name + "_trunc_" + string(truncation_length);
                    end
                    file_name = file_name + "_thick_" + string(double(thick)) + ".mat";
                    load(data_path + file_name, "result");
                    it_repeat = find(result.t_list == median(result.t_list), 1);
                    result_list(it_method).method_name = method;
                    result_list(it_method).info = result.info(it_repeat);
                    if mode == "adaptive"
                        result_list(it_method).method_name = method + ", trunc=" + string(truncation_length);
                    end
                end
            end
            file_name = "qcd_invsqrt_ratios_" + mode + "_sfactor_" + string(sfactor) ...
                + "_thick_" + string(double(thick));
            figure_prefix = figure_path + file_name;
            figure_title = "QCD inverse square root, " + mode + ", sfactor=" + string(sfactor) ...
                + ", thick=" + string(double(thick));
            plot_ratios(result_list, figure_prefix, figure_title, baseline_time);
        end
    end
end
