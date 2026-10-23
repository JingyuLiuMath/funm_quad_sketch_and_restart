exp_qcd_invsqrt_settings;

% Use the same median runs and dimension limit as the QCD results table.
selection.data_path = data_path;
selection.file_prefix = "qcd_invsqrt";
selection.krylov_dim = 150;
selection.thick = false;
group_names = ["Arnoldi", "HmArnoldi"];
method_groups = ["arnoldi", "sarnoldi", "last_orth"; ...
                 "hmarnoldi", "shmarnoldi", "last_hmorth"];
method_labels = ["Arn", "sArn", "adaArn"; ...
                 "HmArn", "sHmArn", "adaHmArn"];

for sfactor = sfactor_list
    selection.sfactor = sfactor;
    for group_index = 1:numel(group_names)
        baseline = read_median_run(selection, "fixed", ...
            method_groups(group_index, 1), inf);
        result_list = struct('method_name', {}, 'group_name', {}, 'info', {});
        for method_index = 1:3
            method = method_groups(group_index, method_index);
            label = method_labels(group_index, method_index);
            if method_index <= 2
                if ~ismember(method, fixed_method_list)
                    continue
                end
                mode = "fixed";
                truncations = inf;
            else
                if selection.krylov_dim ~= adaptive_krylov_dim ...
                        || ~ismember(method, adaptive_update_list)
                    continue
                end
                mode = "adaptive";
                truncations = truncation_length_list;
            end
            for truncation_length = truncations
                info = read_median_run(selection, mode, method, truncation_length);
                method_name = label;
                if mode == "adaptive"
                    method_name = label + "-" + string(truncation_length);
                end
                result_list(end + 1) = struct('method_name', method_name, ...
                    'group_name', group_names(group_index), 'info', info);
            end
        end
        if isempty(result_list)
            continue
        end

        figure_prefix = fullfile(figure_path, "qcd_invsqrt_total_ratios_" ...
            + lower(group_names(group_index)) ...
            + "_krylov_dim_" + string(selection.krylov_dim) ...
            + "_sfactor_" + string(sfactor));
        figure_title = "QCD";
        % Each family uses its own deterministic method as reference.
        plot_ratios_total(result_list, figure_prefix, figure_title, baseline.t_total);
    end
end
