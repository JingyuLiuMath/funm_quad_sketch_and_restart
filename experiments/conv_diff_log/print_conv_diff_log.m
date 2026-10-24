exp_conv_diff_log_settings;

% Select the convection-diffusion results to display.
selection.data_path = data_path;
selection.file_prefix = "conv_diff_log";
selection.thick = false;

% Columns pair the two method families; fixed and adaptive results use separate sections.
method_groups = ["arnoldi", "sarnoldi", "last_orth"; ...
                 "hmarnoldi", "shmarnoldi", "last_hmorth"];
method_labels = ["Arn", "sArn", "adaArn"; ...
                 "HmArn", "sHmArn", "adaHmArn"];
empty_row = struct('method_name', "", 'error', NaN, 'time', NaN, ...
    'n_iter', NaN, 'n_matvec', NaN, 'm_act', NaN);

for sfactor = sfactor_list
    selection.sfactor = sfactor;
    sections = cell(1, numel(krylov_dim_list) + 1);
    section_labels = strings(1, numel(krylov_dim_list) + 1);
    baseline_times = zeros(numel(krylov_dim_list) + 1, 2);

    for dim_index = 1:numel(krylov_dim_list)
        selection.krylov_dim = krylov_dim_list(dim_index);
        rows = repmat(empty_row, 2 + numel(truncation_length_list), 2);

        for group_index = 1:2
            methods = method_groups(group_index, :);
            labels = method_labels(group_index, :);
            baseline = read_median_run(selection, "fixed", methods(1), inf);
            baseline_times(dim_index, group_index) = baseline.t_total;

            % Fixed-dimension methods: Arn/sArn or HmArn/sHmArn.
            for method_index = 1:2
                if ~ismember(methods(method_index), fixed_method_list)
                    continue
                end
                if method_index == 1
                    info = baseline;
                else
                    info = read_median_run(selection, "fixed", methods(method_index), inf);
                end
                rows(method_index, group_index) = ...
                    summarize_run(labels(method_index), info);
            end

            % Adaptive results exist only for the configured maximum dimension.
            if selection.krylov_dim ~= adaptive_krylov_dim ...
                    || ~ismember(methods(3), adaptive_update_list)
                continue
            end
            for trunc_index = 1:numel(truncation_length_list)
                truncation_length = truncation_length_list(trunc_index);
                info = read_median_run(selection, "adaptive", methods(3), truncation_length);
                label = labels(3) + "-" + string(truncation_length);
                rows(2 + trunc_index, group_index) = summarize_run(label, info);
            end
        end

        sections{dim_index} = rows(1:2, :);
        section_labels(dim_index) = "\(m = " + string(selection.krylov_dim) + "\), fixed";
        if selection.krylov_dim == adaptive_krylov_dim
            sections{end} = rows(3:end, :);
            section_labels(end) = "\(m_{\max} = " + string(adaptive_krylov_dim) + "\), adaptive";
            baseline_times(end, :) = baseline_times(dim_index, :);
        end
    end

    table_options.caption = conv_diff_log_caption(sfactor, adaptive_krylov_dim);
    table_options.label = "tab:conv-diff-log-sfactor-" + string(sfactor);
    table_options.output_file = fullfile(figure_path, ...
        "conv_diff_log_table_sfactor_" + string(sfactor) + ".tex");
    table_options.baseline_times = baseline_times;
    table_options.section_labels = section_labels;
    table_options.highlight_fastest = false;
    print_result(sections, table_options);
end

function caption = conv_diff_log_caption(sfactor, adaptive_krylov_dim)
% Keep each caption sentence separate so the wording is easy to edit.
sentences = [ ...
    "Relative error, runtime in seconds, number of restart cycles " + ...
        "\(n_{\mathrm{iter}}\), number of matrix-vector products " + ...
        "\(n_{\mathrm{matvec}}\), and maximum actual Krylov subspace dimension " + ...
        "\(m_{\mathrm{act}}\) attained during each run for the convection-diffusion logarithm example."; ...
    "The column ``speedup'' is the ratio of the runtime of " + ...
        "Arn (left group) or HmArn (right group) to that of the corresponding method, " + ...
        "using the fixed-dimension baseline at \(m=m_{\max}\) for the adaptive section."; ...
    "For fixed-dimension methods, \(m\) is given in the section heading " + ...
        "and \(s=" + string(sfactor) + "m\)."; ...
    "Adaptive methods use \(m_{\max}=" + string(adaptive_krylov_dim) ...
        + "\) and \(\eta=" + string(sfactor) ...
        + "\), and appear in a separate section." ...
];
caption = strjoin(sentences, " ");
end
