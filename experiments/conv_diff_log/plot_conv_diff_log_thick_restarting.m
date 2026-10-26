exp_conv_diff_log_settings;

% Compare independent median-total-time runs, as in the results table.
selection.data_path = data_path;
selection.file_prefix = "conv_diff_log";
selection.krylov_dim = adaptive_krylov_dim;

% Match the error-history figure when both PDFs use .49\textwidth in LaTeX.
plot_options.figure_size = [2300, 1650];
plot_options.axes_position = [480, 400, 1100, 900];
plot_options.font_size = 128;
plot_options.title_font_size = 128;
plot_options.legend_font_size = 104;
plot_options.line_width = 16;
plot_options.marker_size = 16;
plot_options.marker_line_width = 3;
plot_options.legend_token_size = [80, 60];
plot_options.legend_gap = 25;
plot_options.legend_right_padding = 40;
plot_options.title = "Conv-Diff";
% Focus on the first 10000 matvecs, retaining every run's complete history.
plot_options.x_limits = [0, 10000];
% Keep the full-count tick labels separate at the reference font size.
plot_options.x_ticks = 0:5000:10000;
% Match the graph-Laplacian thick-restarting figure's error scale.
plot_options.y_limits = [1e-12, 1];
plot_options.y_ticks = 10 .^ (-12:2:0);

for sfactor = sfactor_list
    selection.sfactor = sfactor;
    info_without = {};
    info_with = {};
    method_names = strings(1, 0);
    if selection.krylov_dim == adaptive_krylov_dim ...
            && ismember("last_orth", adaptive_update_list)
        for truncation_length = truncation_length_list
            [info_without{end + 1}, info_with{end + 1}] = ...
                read_run_pair(selection, "adaptive", "last_orth", truncation_length); %#ok<SAGROW>
            method_names(end + 1) = "adaArn-" + string(truncation_length); %#ok<SAGROW>
        end
    end
    if isempty(method_names)
        continue
    end

    figure_name = "conv_diff_log_thick_restarting_arnoldi" ...
        + "_krylov_dim_" + string(selection.krylov_dim) ...
        + "_sfactor_" + string(sfactor);
    figure_prefix = fullfile(figure_path, figure_name);
    plot_conv_thick_restarting(info_without, info_with, figure_prefix, ...
        method_names, plot_options);
end

function [info_without, info_with] = read_run_pair(selection, mode, method, truncation_length)
selection.thick = false;
info_without = read_conv_thick_run(selection, mode, method, truncation_length);
selection.thick = true;
info_with = read_conv_thick_run(selection, mode, method, truncation_length);
end

function info = read_conv_thick_run(selection, mode, method, truncation_length)
% Keep this figure independent of shared plotting and result-loading helpers.
factor = selection.sfactor;
if mode == "fixed" && ismember(method, ["arnoldi", "hmarnoldi"])
    factor = 0;
end
file_name = selection.file_prefix + "_" + mode + "_" + method ...
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

function fig = plot_conv_thick_restarting(info_without, info_with, figure_prefix, method_names, options)
% Every line retains its complete history; xlim hides long non-thick tails.
if nargin < 5
    options = struct();
end
method_names = string(method_names);
method_names = method_names(:);
num_methods = numel(method_names);
assert(iscell(info_without) && iscell(info_with) && num_methods >= 1 && ...
    numel(info_without) == num_methods && numel(info_with) == num_methods, ...
    'Supply matching info cell arrays and adaptive method names.');
infos = [info_without(:), info_with(:)];
font_size = conv_thick_option(options, 'font_size', 128);
legend_font_size = conv_thick_option(options, 'legend_font_size', 104);
title_font_size = conv_thick_option(options, 'title_font_size', font_size);
line_width = conv_thick_option(options, 'line_width', 16);
marker_size = conv_thick_option(options, 'marker_size', 16);
marker_line_width = conv_thick_option(options, 'marker_line_width', 3);
figure_size = conv_thick_option(options, 'figure_size', [2300, 1650]);
axes_position = conv_thick_option(options, 'axes_position', [480, 400, 1100, 900]);
% Use the same adaptive colors and marker shapes as the error-history figure.
colors = lines(num_methods + 2);
colors = colors(3:end, :);
markers = ["^", "d", "v", "o", "s", ">", "<", "p"];
styles = {'-', '--'};
x_histories = cell(num_methods, 2);
error_histories = cell(num_methods, 2);
positive_errors = [];
for method_index = 1:num_methods
    for restart_index = 1:2
        info = infos{method_index, restart_index};
        assert(isstruct(info) && isscalar(info) && ...
            all(isfield(info, {'num_oracle', 'err_history'})), ...
            'Each run must contain num_oracle and err_history.');
        counts = info.num_oracle(:);
        errors = info.err_history(:);
        assert(isnumeric(counts) && isreal(counts) && ~isempty(counts) && ...
            all(isfinite(counts)) && all(counts >= 0) && ...
            isnumeric(errors) && isreal(errors) && numel(counts) == numel(errors), ...
            'Oracle counts and error history must be nonempty matching vectors.');
        x_histories{method_index, restart_index} = cumsum(counts);
        positive_errors = [positive_errors; ...
            errors(isfinite(errors) & errors > 0)]; %#ok<AGROW>
        errors(~isfinite(errors) | errors <= 0) = NaN;
        error_histories{method_index, restart_index} = errors;
    end
end
assert(~isempty(positive_errors), 'At least one error must be positive and finite.');
y_limits = conv_thick_option(options, 'y_limits', [1e-12, 1]);
y_ticks = conv_thick_option(options, 'y_ticks', 10 .^ (-12:2:0));
x_limits = conv_thick_option(options, 'x_limits', [0, 10000]);
x_ticks = conv_thick_option(options, 'x_ticks', 0:5000:10000);

fig = figure('Color', 'white', 'Units', 'pixels', ...
    'Position', [100, 100, figure_size], 'Renderer', 'painters');
ax = axes(fig, 'Units', 'pixels', 'Position', axes_position, ...
    'PositionConstraint', 'innerposition', 'FontSize', font_size, ...
    'FontWeight', 'normal', 'TickLabelInterpreter', 'latex', ...
    'Box', 'on', 'TickDir', 'in', 'Layer', 'top', ...
    'XMinorTick', 'off', 'YMinorTick', 'on', 'YScale', 'log', ...
    'YLim', y_limits, 'YTick', y_ticks, 'Tag', 'thick_focused');
hold(ax, 'on');
line_handles = gobjects(num_methods, 2);
for method_index = 1:num_methods
    for restart_index = 1:2
        line_handles(method_index, restart_index) = semilogy(ax, ...
            x_histories{method_index, restart_index}, ...
            error_histories{method_index, restart_index}, ...
            'Color', colors(method_index, :), 'LineStyle', styles{restart_index}, ...
            'LineWidth', line_width, 'Marker', 'none', ...
            'HandleVisibility', 'off', 'Tag', 'restart_curve');
    end
end
% Use about twelve markers per visible curve, as in the error-history plot.
% Offset the methods' marker positions to reduce overlap.
legend_handles = gobjects(num_methods, 1);
for method_index = 1:num_methods
    marker = markers(1 + mod(method_index - 1, numel(markers)));
    for restart_index = 1:2
        visible_cycles = find(x_histories{method_index, restart_index} >= x_limits(1) ...
            & x_histories{method_index, restart_index} <= x_limits(2));
        marker_step = max(1, ceil(numel(visible_cycles) / 12));
        marker_start = 1 + floor((method_index - 1) * marker_step / num_methods);
        marker_indices = visible_cycles(marker_start:marker_step:end);
        semilogy(ax, x_histories{method_index, restart_index}, ...
            error_histories{method_index, restart_index}, ...
            'Color', colors(method_index, :), 'LineStyle', 'none', ...
            'LineWidth', marker_line_width, 'Marker', marker, ...
            'MarkerSize', marker_size, 'MarkerFaceColor', 'white', ...
            'MarkerIndices', marker_indices, ...
            'HandleVisibility', 'off', 'Tag', 'cycle_markers');
    end
    % Keep method symbols legible in the legend, as in the reference figure.
    legend_handles(method_index) = semilogy(ax, NaN, NaN, ...
        'Color', colors(method_index, :), 'LineStyle', '-', ...
        'LineWidth', 12, 'Marker', marker, 'MarkerSize', 40, ...
        'MarkerFaceColor', 'none', 'Tag', 'method_legend');
end
set(ax, 'YScale', 'log', 'YLim', y_limits, 'YTick', y_ticks);
xlim(ax, x_limits);
xticks(ax, x_ticks);
xtickangle(ax, 0);
% Show full counts so MATLAB's scientific-notation multiplier cannot cover the xlabel.
ax.XAxis.Exponent = 0;
grid(ax, 'off');
xlabel(ax, 'number of matrix-vector products', 'Interpreter', 'latex', ...
    'FontSize', font_size, 'FontWeight', 'normal');
ylabel(ax, 'relative error', 'Interpreter', 'latex', ...
    'FontSize', font_size, 'FontWeight', 'normal');
title_text = conv_thick_option(options, 'title', "Conv-Diff");
title_text = replace(replace(string(title_text), "_", "\_"), " ", "\ ");
title(ax, "$\mathrm{" + title_text + "}$", 'Interpreter', 'latex', ...
    'FontSize', title_font_size, 'FontWeight', 'normal');
lgd = legend(ax, legend_handles, replace(method_names, "_", "\_"), ...
    'Interpreter', 'latex', 'FontSize', legend_font_size, 'FontWeight', 'normal', ...
    'NumColumns', 1, 'Box', 'on', 'Location', 'eastoutside', 'AutoUpdate', 'off');
lgd.ItemTokenSize = conv_thick_option(options, 'legend_token_size', [80, 60]);
lgd.Units = 'pixels';

% Preserve the same plotting area and label margins as the error-history PDF.
drawnow;
ax.Position = axes_position;
drawnow;
inset = ax.TightInset;
axes_position(1:2) = max(axes_position(1:2), inset(1:2) + 25);
ax.Position = axes_position;
fig.Position(3:4) = max(fig.Position(3:4), ...
    axes_position(1:2) + axes_position(3:4) + inset(3:4) + 25);
drawnow;

% An outside legend must not shrink the requested plotting area.
ax.Position = axes_position;
legend_position = lgd.Position;
legend_position(1) = axes_position(1) + axes_position(3) ...
    + conv_thick_option(options, 'legend_gap', 25);
legend_position(2) = axes_position(2) + (axes_position(4) - legend_position(4)) / 2;
lgd.Position = legend_position;
fig.Position(3) = max(fig.Position(3), legend_position(1) + legend_position(3) ...
    + conv_thick_option(options, 'legend_right_padding', 40));
drawnow;

figure_prefix = string(figure_prefix);
figure_directory = fileparts(figure_prefix);
if strlength(figure_directory) > 0 && ~isfolder(figure_directory)
    mkdir(figure_directory);
end
paper_size = fig.Position(3:4) / get(groot, 'ScreenPixelsPerInch');
set(fig, 'PaperUnits', 'inches', 'PaperSize', paper_size, ...
    'PaperPosition', [0, 0, paper_size], 'PaperPositionMode', 'manual');
print(fig, char(figure_prefix + ".pdf"), '-dpdf', '-painters');
print(fig, char(figure_prefix + ".png"), '-dpng', '-r300');
end

function value = conv_thick_option(options, name, default_value)
if isfield(options, name)
    value = options.(name);
else
    value = default_value;
end
end
