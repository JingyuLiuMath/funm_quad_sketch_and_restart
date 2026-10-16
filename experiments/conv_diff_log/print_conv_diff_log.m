exp_conv_diff_log_settings;

mode_list = ["fixed", "adaptive"];
for sfactor = sfactor_list
    for truncation_length = truncation_length_list
        fprintf('\nsfactor = %d, trunc = %d (fixed results reused; unsketched sfactor=0)\n', ...
            sfactor, truncation_length);
        fprintf('%-9s %-14s %7s %7s %5s %12s %12s %6s %9s\n', ...
            'mode', 'method', 'sfactor', 'trunc', 'thick', 'time', 'true error', 'iter', 'oracle');
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
                    file_name = "conv_diff_log_" + mode + "_" + method ...
                        + "_sfactor_" + string(factor);
                    if mode == "adaptive"
                        file_name = file_name + "_trunc_" + string(truncation_length);
                    end
                    file_name = file_name + "_thick_" + string(double(thick)) + ".mat";
                    load(data_path + file_name, "result");
                    it_repeat = find(result.t_list == median(result.t_list), 1);
                    info = result.info(it_repeat);
                    fprintf('%-9s %-14s %7d %7g %5d %12.5e %12.5e %6d %9d\n', ...
                        mode, method, result.sfactor, result.param.truncation_length, thick, ...
                        result.t_list(it_repeat), result.err_list(it_repeat), ...
                        info.iter, sum(info.num_oracle));
                end
            end
        end
    end
end
