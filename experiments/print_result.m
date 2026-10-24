function latex_text = print_result(rows, options)
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
            cells = format_result(result, options.baseline_times(section_index, group));
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

function cells = format_result(result, baseline_time)
speedup = baseline_time / result.time;
cells = [ ...
    string(result.method_name), ...
    latex_scientific(result.error), ...
    string(sprintf('%.1f', result.time)), ...
    string(result.n_iter), ...
    string(result.n_matvec), ...
    string(result.m_act), ...
    "\(" + string(sprintf('%.1f', speedup)) + "\,\times\)" ...
];
end

function text = latex_scientific(value)
if ~isfinite(value)
    text = "\(\mathrm{" + string(value) + "}\)";
    return
end
parts = split(string(sprintf('%.1e', value)), 'e');
text = parts(1) + "e" + string(str2double(parts(2)));
end
