exp_qcd_invsqrt_settings;

% Use the same selected runs and method names as the QCD table.
selection.data_path = data_path;
selection.file_prefix = "qcd_invsqrt";
selection.krylov_dim = 150;
selection.thick = false;
group_names = ["arnoldi", "hmarnoldi"];
method_groups = ["arnoldi", "sarnoldi", "last_orth"; ...
                 "hmarnoldi", "shmarnoldi", "last_hmorth"];
method_labels = ["Arn", "sArn", "adaArn"; ...
                 "HmArn", "sHmArn", "adaHmArn"];

% Enlarge the plotting area for two figures placed side by side in LaTeX.
plot_options.figure_size = [2300, 1650];
plot_options.axes_position = [480, 400, 1100, 900];
plot_options.font_size = 128;
plot_options.legend_font_size = 104;
plot_options.line_width = 12;
plot_options.marker_size = 40;
plot_options.legend_location = "eastoutside";
plot_options.legend_gap = 25;
plot_options.legend_token_size = [80, 60];
plot_options.x_label = "number of matrix-vector products";
plot_options.y_label = "relative error";
plot_options.x_limits = [0, 1500];
plot_options.x_ticks = 0:500:1500;
% Leave room below the last markers while retaining the decade tick labels.
plot_options.y_limits = [3e-13, 1];

for sfactor = sfactor_list
    selection.sfactor = sfactor;
    for group_index = 1:numel(group_names)
        result_list = struct('method_name', {}, 'res_hist', {}, 'matvec_hist', {});
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
                result.method_name = label;
                if mode == "adaptive"
                    result.method_name = label + "-" + string(truncation_length);
                end
                result.res_hist = info.err_history(:).';
                result.matvec_hist = cumsum(info.num_oracle(:)).';
                result_list(end + 1) = result;
            end
        end
        if isempty(result_list)
            continue
        end
        figure_prefix = fullfile(figure_path, "qcd_invsqrt_err_" ...
            + group_names(group_index) + "_sfactor_" + string(sfactor));
        plot_res(result_list, figure_prefix, "QCD", plot_options);
    end
end
