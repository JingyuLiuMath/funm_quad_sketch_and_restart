exp_graph_laplacian_invsqrt_settings;

% Use the same configurations and median runs as the thick-restart plots.
selection.data_path = data_path;
selection.file_prefix = "frac_laplacian_invsqrt";
fixed_methods = ["arnoldi", "sarnoldi"];
fixed_labels = ["Arn", "sArn"];
adaptive_methods = ["last_orth", "last_hmorth"];
adaptive_labels = ["adaArn", "adaHmArn"];

% Use the same m=150 selection as the comparison plots.
for krylov_dim = 150
    selection.krylov_dim = krylov_dim;
    for sfactor = sfactor_list
        selection.sfactor = sfactor;
        rows = struct('method_name', {}, 'n_iter', {}, ...
            'n_matvec', {}, 'speedup', {});

        % Fixed harmonic methods do not support thick restarting.
        for method_index = 1:numel(fixed_methods)
            method = fixed_methods(method_index);
            if ~ismember(method, fixed_method_list)
                continue
            end
            rows(end + 1) = thick_restart_row(selection, "fixed", method, ...
                fixed_labels(method_index), inf); %#ok<SAGROW>
        end

        if krylov_dim == adaptive_krylov_dim
            for method_index = 1:numel(adaptive_methods)
                method = adaptive_methods(method_index);
                if ~ismember(method, adaptive_update_list)
                    continue
                end
                for truncation_length = truncation_length_list
                    label = adaptive_labels(method_index) + "-" ...
                        + string(truncation_length);
                    rows(end + 1) = thick_restart_row(selection, "adaptive", ...
                        method, label, truncation_length); %#ok<SAGROW>
                end
            end
        end
        if isempty(rows)
            continue
        end

        output_file = fullfile(figure_path, ...
            "graph_laplacian_invsqrt_thick_restarting_table_krylov_dim_" ...
            + string(krylov_dim) + "_sfactor_" + string(sfactor) + ".tex");
        print_thick_restart_table(rows, selection, output_file);
    end
end

function row = thick_restart_row(selection, mode, method, label, truncation_length)
selection.thick = false;
info_without = graph_read_median_run(selection, mode, method, truncation_length);
selection.thick = true;
info_with = graph_read_median_run(selection, mode, method, truncation_length);

row = struct('method_name', label, 'n_iter', info_with.iter, ...
    'n_matvec', sum(info_with.num_oracle(:)), ...
    'speedup', info_without.t_total / info_with.t_total);
end

function print_thick_restart_table(rows, selection, output_file)
caption = "Restart cycles \(n_{\mathrm{iter}}\) and matrix-vector products " ...
    + "\(n_{\mathrm{matvec}}\) with thick restarting for the graph Laplacian " ...
    + "example. Speedup is the total runtime without thick restarting divided " ...
    + "by that with thick restarting for the same method and configuration. " ...
    + "Each run is selected by its median total runtime. The fixed Krylov " ...
    + "dimension or adaptive dimension limit is " + string(selection.krylov_dim) ...
    + ", and the sketching size factor is " + string(selection.sfactor) + ".";
table_label = "tab:graph-laplacian-invsqrt-thick-restarting-krylov-dim-" ...
    + string(selection.krylov_dim) + "-sfactor-" + string(selection.sfactor);
latex_lines = [ ...
    "% Requires \usepackage{booktabs} in the preamble."; ...
    "\begin{table}[htbp]"; ...
    "\centering"; ...
    "\caption{" + caption + "}"; ...
    "\label{" + table_label + "}"; ...
    "\begin{tabular}{@{}lrrr@{}}"; ...
    "\toprule"; ...
    "method & \(n_{\mathrm{iter}}\) & \(n_{\mathrm{matvec}}\) & speedup \\"; ...
    "\midrule" ...
];
data_lines = strings(numel(rows), 1);
for row_index = 1:numel(rows)
    row = rows(row_index);
    cells = [replace(string(row.method_name), "_", "\_"), ...
        string(row.n_iter), string(row.n_matvec), ...
        "\(" + string(sprintf('%.1f', row.speedup)) + "\,\times\)"];
    data_lines(row_index) = strjoin(cells, " & ") + " \\";
end
latex_lines = [latex_lines; data_lines; "\bottomrule"; "\end{tabular}"; "\end{table}"];
latex_text = strjoin(latex_lines, newline) + newline;
fprintf('%s', latex_text);

fid = fopen(output_file, 'w');
assert(fid ~= -1, 'Cannot open LaTeX table output file: %s', output_file);
file_cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s', latex_text);
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
