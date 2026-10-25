exp_graph_laplacian_invsqrt_settings;

% Select the graph Laplacian results to display.
selection.data_path = data_path;
selection.file_prefix = "frac_laplacian_invsqrt";
selection.krylov_dim = 150;
selection.thick = false;

% Columns correspond to the two method families; rows pair matching methods.
method_groups = ["arnoldi", "sarnoldi", "last_orth"; ...
                 "hmarnoldi", "shmarnoldi", "last_hmorth"];
method_labels = ["Arn", "sArn", "adaArn"; ...
                 "HmArn", "sHmArn", "adaHmArn"];
empty_row = struct('method_name', "", 'error', NaN, 'time', NaN, ...
    'n_iter', NaN, 'n_matvec', NaN, 'm_act', NaN);

for sfactor = sfactor_list
    selection.sfactor = sfactor;
    rows = repmat(empty_row, 2 + numel(truncation_length_list), 2);
    baseline_times = zeros(1, 2);

    for group_index = 1:2
        methods = method_groups(group_index, :);
        labels = method_labels(group_index, :);
        baseline = read_median_run(selection, "fixed", methods(1), inf);
        baseline_times(group_index) = baseline.t_total;

        % Fixed-dimension methods: Arn/sArn or HmArn/sHmArn.
        for method_index = 1:2
            if ~ismember(methods(method_index), fixed_method_list)
                continue
            end
            if method_index == 1
                info = baseline;
            elseif method_index == 2
                info = read_median_run(selection, "fixed", methods(method_index), inf);
            end
            rows(method_index, group_index) = ...
                summarize_run(labels(method_index), info);
        end

        % Adaptive methods: one row for each truncation length.
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

    % The shared printer handles formatting, speedup, highlighting and output.
    table_options.caption = graph_laplacian_caption(selection);
    table_options.label = "tab:graph-laplacian-invsqrt-sfactor-" + string(sfactor);
    table_options.output_file = fullfile(figure_path, ...
        "graph_laplacian_invsqrt_table_sfactor_" + string(sfactor) + ".tex");
    table_options.baseline_times = baseline_times;
    print_result(rows, table_options);
end

function caption = graph_laplacian_caption(selection)
% Keep each caption sentence separate so the wording is easy to edit.
sentences = [ ...
    "Relative error, runtime in seconds, number of restart cycles " + ...
        "\(n_{\mathrm{iter}}\), number of matrix-vector products " + ...
        "\(n_{\mathrm{matvec}}\), and average actual Krylov subspace dimension " + ...
        "\(m_{\mathrm{act}}\) per restart cycle for the graph Laplacian example."; ...
    "The column ``speedup'' is defined as the ratio of the runtime of " + ...
        "Arn (left group) or HmArn (right group) to that of the corresponding method."; ...
    "For fixed-dimension methods, \(m=" + string(selection.krylov_dim) ...
        + "\) and \(s=" + string(selection.sfactor) + "m=" ...
        + string(selection.sfactor * selection.krylov_dim) + "\)."; ...
    "For adaptive methods, \(m_{\max}=" + string(selection.krylov_dim) ...
        + "\) and \(\eta=" + string(selection.sfactor) + "\)."; ...
    "The row corresponding to the fastest method in each group is highlighted in bold." ...
];
caption = strjoin(sentences, " ");
end
