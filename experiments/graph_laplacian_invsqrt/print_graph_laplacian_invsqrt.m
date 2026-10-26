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
        baseline = graph_read_median_run(selection, "fixed", methods(1), inf);
        baseline_times(group_index) = baseline.t_total;

        % Fixed-dimension methods: Arn/sArn or HmArn/sHmArn.
        for method_index = 1:2
            if ~ismember(methods(method_index), fixed_method_list)
                continue
            end
            if method_index == 1
                info = baseline;
            elseif method_index == 2
                info = graph_read_median_run(selection, "fixed", methods(method_index), inf);
            end
            rows(method_index, group_index) = ...
                graph_summarize_run(labels(method_index), info);
        end

        % Adaptive methods: one row for each truncation length.
        if selection.krylov_dim ~= adaptive_krylov_dim ...
                || ~ismember(methods(3), adaptive_update_list)
            continue
        end
        for trunc_index = 1:numel(truncation_length_list)
            truncation_length = truncation_length_list(trunc_index);
            info = graph_read_median_run(selection, "adaptive", methods(3), truncation_length);
            label = labels(3) + "-" + string(truncation_length);
            rows(2 + trunc_index, group_index) = graph_summarize_run(label, info);
        end
    end

    % The local printer handles formatting, speedup, highlighting and output.
    table_options.caption = graph_laplacian_caption(selection);
    table_options.label = "tab:graph-laplacian-invsqrt-sfactor-" + string(sfactor);
    table_options.output_file = fullfile(figure_path, ...
        "graph_laplacian_invsqrt_table_sfactor_" + string(sfactor) + ".tex");
    table_options.baseline_times = baseline_times;
    graph_print_result(rows, table_options);
end

function caption = graph_laplacian_caption(selection)
% Keep each caption sentence separate so the wording is easy to edit.
sentences = [ ...
    "Relative error, runtime in seconds, number of restart cycles " + ...
        "\(n_{\mathrm{iter}}\), number of matrix-vector products " + ...
        "\(n_{\mathrm{matvec}}\), and average actual Krylov subspace dimension " + ...
        "\(m_{\mathrm{act}}\) per restart cycle for the graph Laplacian example."; ...
    "The column ``speedup'' reports the percentage runtime reduction " + ...
        "\((t_{\mathrm{base}}-t_{\mathrm{curr}})/t_{\mathrm{base}}\times100\%\), " + ...
        "where the baseline is Arn (left group) or HmArn (right group)."; ...
    "For fixed-dimension methods, \(m=" + string(selection.krylov_dim) ...
        + "\) and \(s=" + string(selection.sfactor) + "m=" ...
        + string(selection.sfactor * selection.krylov_dim) + "\)."; ...
    "For adaptive methods, \(m_{\max}=" + string(selection.krylov_dim) ...
        + "\) and \(\eta=" + string(selection.sfactor) + "\)."; ...
    "The row corresponding to the fastest method in each group is highlighted in bold." ...
];
caption = strjoin(sentences, " ");
end

% Graph-local result loading and table formatting.

function info = graph_read_median_run(selection, mode, method, truncation_length)
%READ_MEDIAN_RUN Load the run with median total runtime for a configuration.
% selection supplies data_path, file_prefix, krylov_dim, sfactor, and thick.
% All reported values come from the same selected repeat.

mode = string(mode);
method = string(method);
factor = selection.sfactor;
if mode == "fixed" && ismember(method, ["arnoldi", "hmarnoldi"])
    factor = 0;
end
file_name = string(selection.file_prefix) + "_" + mode + "_" + method ...
    + "_krylov_dim_" + string(selection.krylov_dim) ...
    + "_sfactor_" + string(factor) ...
    + "_thick_" + string(double(selection.thick));
if mode == "adaptive"
    file_name = file_name + "_trunc_" + string(truncation_length);
end
data = load(fullfile(selection.data_path, file_name + ".mat"), "result");
repeat_index = find(data.result.t_list == median(data.result.t_list), 1);
info = data.result.info(repeat_index);
end

function row = graph_summarize_run(label, info)
%SUMMARIZE_RUN Normalize a selected run for the shared results table.
% m_act is the mean actual Krylov dimension across restart cycles.

row = struct('method_name', string(label), 'error', info.final_err, ...
    'time', info.t_total, 'n_iter', info.iter, ...
    'n_matvec', sum(info.num_oracle(:)), ...
    'm_act', mean(info.dim(:)));
end

function latex_text = graph_print_result(rows, options)
%PRINT_RESULT Print and save a LaTeX table of grouped experiment results.
% rows is an N-by-G struct array, with one method per row in each group.
% Each entry has method_name, error, time, n_iter, n_matvec, and m_act;
% an empty method_name denotes a missing result. options supplies caption,
% label, output_file, and the unrounded baseline_times for the G groups.
% For multiple sections, rows is a cell array of N-by-G struct arrays,
% options.section_labels names each section, and baseline_times is S-by-G.
% Set options.highlight_fastest to false to disable bolding the fastest methods.

if ~iscell(rows)
    rows = {rows};
end
num_groups = size(rows{1}, 2);
highlight_fastest = ~isfield(options, 'highlight_fastest') || options.highlight_fastest;

group_header = "method & error & time & \(n_{\mathrm{iter}}\) & " + ...
    "\(n_{\mathrm{matvec}}\) & \(m_{\mathrm{act}}\) & speedup";
column_spec = "@{}" + strjoin(repmat("lrrrrrr", 1, num_groups), ...
    "@{\quad}|@{\quad}") + "@{}";
latex_lines = [ ...
    "% Requires \usepackage{graphicx,booktabs} in the preamble."; ...
    "\begin{table}[htbp]"; ...
    "\centering"; ...
    "\caption{" + string(options.caption) + "}"; ...
    "\label{" + string(options.label) + "}"; ...
    "\resizebox{\textwidth}{!}{%"; ...
    "\begin{tabular}{" + column_spec + "}"; ...
    "\toprule"; ...
    strjoin(repmat(group_header, 1, num_groups), " & ") + " \\"; ...
    "\midrule" ...
];

for section_index = 1:numel(rows)
    section_rows = rows{section_index};
    has_result = arrayfun(@(row) strlength(string(row.method_name)) > 0, section_rows);
    keep_rows = any(has_result, 2);
    section_rows = section_rows(keep_rows, :);
    has_result = has_result(keep_rows, :);

    % Compare original runtimes within this section, before display rounding.
    fastest_times = nan(1, num_groups);
    for group = 1:num_groups
        times = [section_rows(has_result(:, group), group).time];
        times = times(isfinite(times));
        if ~isempty(times)
            fastest_times(group) = min(times);
        end
    end

    if section_index > 1
        latex_lines(end + 1, 1) = "\midrule";
    end
    if isfield(options, 'section_labels')
        latex_lines(end + 1, 1) = "\multicolumn{" + string(7 * num_groups) ...
            + "}{c}{" + string(options.section_labels(section_index)) + "} \\";
        latex_lines(end + 1, 1) = "\midrule";
    end

    for row_index = 1:size(section_rows, 1)
        row_cells = repmat("--", 1, 7 * num_groups);
        for group = 1:num_groups
            if ~has_result(row_index, group)
                continue
            end
            result = section_rows(row_index, group);
            cells = graph_format_result(result, options.baseline_times(section_index, group));
            if highlight_fastest && result.time == fastest_times(group)
                % A group-local winner has every cell in bold, including math.
                cells = "{\bfseries\boldmath " + cells + "}";
            end
            row_cells(7 * (group - 1) + (1:7)) = cells;
        end
        latex_lines(end + 1, 1) = strjoin(row_cells, " & ") + " \\";
    end
end

latex_lines = [latex_lines; "\bottomrule"; "\end{tabular}%"; ...
    "}"; "\end{table}"];
latex_text = strjoin(latex_lines, newline) + newline;
fprintf('%s', latex_text);

fid = fopen(options.output_file, 'w');
assert(fid ~= -1, 'Cannot open LaTeX table output file: %s', options.output_file);
file_cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '%s', latex_text);
end

function cells = graph_format_result(result, baseline_time)
speedup = (baseline_time - result.time) / baseline_time;
cells = [ ...
    string(result.method_name), ...
    graph_latex_scientific(result.error), ...
    string(sprintf('%.1f', result.time)), ...
    string(result.n_iter), ...
    string(result.n_matvec), ...
    string(sprintf('%.1f', result.m_act)), ...
    "\(" + string(sprintf('%.1f', 100 * speedup)) + "\%\)" ...
];
end

function text = graph_latex_scientific(value)
if ~isfinite(value)
    text = "\(\mathrm{" + string(value) + "}\)";
    return
end
parts = split(string(sprintf('%.1e', value)), 'e');
text = parts(1) + "e" + string(str2double(parts(2)));
end
