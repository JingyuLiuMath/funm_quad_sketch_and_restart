exp_graph_laplacian_invsqrt_settings;

% Compare independent median-total-time runs, as in the results table.
selection.data_path = data_path;
selection.file_prefix = "frac_laplacian_invsqrt";
selection.krylov_dim = 150;
fixed_methods = ["arnoldi", "sarnoldi"];
fixed_labels = ["Arn", "sArn"];

% Focus on the range containing every thick-restart and fixed-method run.
plot_options.figure_size = [3200, 2450];
plot_options.font_size = 128;
plot_options.title_font_size = 128;
plot_options.legend_font_size = 104;
plot_options.line_width = 12;
plot_options.legend_token_size = [80, 60];
plot_options.legend_num_columns = 2;
plot_options.title = "Graph Laplacian";
plot_options.subtitle = "solid: without thick; dashed: with thick";
plot_options.subtitle_font_size = 104;
% Optional ranges for the view and the error axis:
% plot_options.zoom_x_limits = [0, 1000];
% plot_options.y_limits = [1e-8, 1e-4];

for sfactor = sfactor_list
    selection.sfactor = sfactor;
    info_without = {};
    info_with = {};
    method_names = strings(1, 0);
    for method_index = 1:numel(fixed_methods)
        method = fixed_methods(method_index);
        if ~ismember(method, fixed_method_list)
            continue
        end
        [info_without{end + 1}, info_with{end + 1}] = ...
            read_run_pair(selection, "fixed", method, inf); %#ok<SAGROW>
        method_names(end + 1) = fixed_labels(method_index); %#ok<SAGROW>
    end

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

    figure_name = "graph_laplacian_invsqrt_thick_restarting_arnoldi" ...
        + "_krylov_dim_" + string(selection.krylov_dim) ...
        + "_sfactor_" + string(sfactor);
    figure_prefix = fullfile(figure_path, figure_name);
    plot_graph_thick_focused(info_without, info_with, figure_prefix, ...
        method_names, plot_options);
end

function [info_without, info_with] = read_run_pair(selection, mode, method, truncation_length)
selection.thick = false;
info_without = read_graph_thick_run(selection, mode, method, truncation_length);
selection.thick = true;
info_with = read_graph_thick_run(selection, mode, method, truncation_length);
end

function info = read_graph_thick_run(selection, mode, method, truncation_length)
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

function fig = plot_graph_thick_focused(info_without, info_with, figure_prefix, method_names, options)
% The visible range covers all thick runs and both fixed non-thick runs.
% Every line retains its complete history; xlim hides long adaptive tails.
if nargin < 5
    options = struct();
end
method_names = string(method_names);
method_names = method_names(:);
num_methods = numel(method_names);
assert(iscell(info_without) && iscell(info_with) && num_methods >= 2 && ...
    numel(info_without) == num_methods && numel(info_with) == num_methods, ...
    'Supply matching info cell arrays and method names, with fixed methods first.');
infos = [info_without(:), info_with(:)];
font_size = graph_thick_option(options, 'font_size', 128);
legend_font_size = graph_thick_option(options, 'legend_font_size', 104);
title_font_size = graph_thick_option(options, 'title_font_size', font_size);
subtitle_font_size = graph_thick_option(options, 'subtitle_font_size', legend_font_size);
line_width = graph_thick_option(options, 'line_width', 12);
panel_size = graph_thick_option(options, 'panel_size', [2200, 1000]);
figure_size = max([3200, 2400], ...
    graph_thick_option(options, 'figure_size', [3200, 2400]));
colors = lines(num_methods);
styles = {'-', '--'};
x_histories = cell(num_methods, 2);
error_histories = cell(num_methods, 2);
legend_labels = strings(1, num_methods);
positive_errors = [];
zoom_xmax = 0;
for method_index = 1:num_methods
    for restart_index = 1:2
        info = infos{method_index, restart_index};
        assert(isstruct(info) && isscalar(info) && ...
            all(isfield(info, {'num_oracle', 'err_history', 't_total'})), ...
            'Each run must contain num_oracle, err_history, and t_total.');
        counts = info.num_oracle(:);
        errors = info.err_history(:);
        assert(isnumeric(counts) && isreal(counts) && ~isempty(counts) && ...
            all(isfinite(counts)) && all(counts >= 0) && ...
            isnumeric(errors) && isreal(errors) && numel(counts) == numel(errors), ...
            'Oracle counts and error history must be nonempty matching vectors.');
        assert(isnumeric(info.t_total) && isscalar(info.t_total) && ...
            isreal(info.t_total) && isfinite(info.t_total) && info.t_total >= 0);
        x_histories{method_index, restart_index} = cumsum(counts);
        positive_errors = [positive_errors; ...
            errors(isfinite(errors) & errors > 0)]; %#ok<AGROW>
        errors(~isfinite(errors) | errors <= 0) = NaN;
        error_histories{method_index, restart_index} = errors;
        if restart_index == 2 || method_index <= 2
            zoom_xmax = max(zoom_xmax, sum(counts));
        end
    end
    legend_labels(method_index) = sprintf('%s (%.1f / %.1f s)', ...
        char(method_names(method_index)), infos{method_index, 1}.t_total, ...
        infos{method_index, 2}.t_total);
end
assert(~isempty(positive_errors), 'At least one error must be positive and finite.');
y_limits = graph_thick_option(options, 'y_limits', ...
    [min(positive_errors) / 3, max(1, 1.1 * max(positive_errors))]);
exponents = [ceil(log10(y_limits(1))), floor(log10(y_limits(2)))];
exponent_step = max(1, ceil(diff(exponents) / 6));
y_ticks = graph_thick_option(options, 'y_ticks', ...
    10 .^ (exponents(1):exponent_step:exponents(2)));
[zoom_limits, zoom_ticks] = graph_thick_zoom_range(zoom_xmax);
zoom_limits = graph_thick_option(options, 'zoom_x_limits', zoom_limits);
zoom_ticks = graph_thick_option(options, 'zoom_x_ticks', zoom_ticks);

fig = figure('Color', 'white', 'Units', 'pixels', ...
    'Position', [100, 100, figure_size], 'Renderer', 'painters');
axes_position = [500, 800, panel_size];
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
            'LineWidth', line_width, 'Marker', 'none');
    end
end
set(ax, 'YScale', 'log', 'YLim', y_limits, 'YTick', y_ticks);
xlim(ax, zoom_limits);
xticks(ax, zoom_ticks);
xtickangle(ax, 0);
grid(ax, 'off');
ylabel(ax, 'relative error', 'Interpreter', 'latex', ...
    'FontSize', font_size, 'FontWeight', 'normal');
lgd = legend(ax, line_handles(:, 1), replace(legend_labels, "_", "\_"), ...
    'Interpreter', 'latex', 'FontSize', legend_font_size, 'FontWeight', 'normal', ...
    'NumColumns', graph_thick_option(options, 'legend_num_columns', 2), ...
    'Box', 'off', 'Location', 'southoutside', 'AutoUpdate', 'off');
lgd.ItemTokenSize = graph_thick_option(options, 'legend_token_size', [100, 60]);
lgd.Units = 'pixels';

canvas_ax = axes(fig, 'Units', 'pixels', 'Position', [0, 0, figure_size], ...
    'PositionConstraint', 'innerposition', 'Color', 'none', 'Visible', 'off', ...
    'HandleVisibility', 'off', 'Tag', 'thick_canvas_labels');
title_text = graph_thick_option(options, 'title', "Graph Laplacian");
title_text = replace(replace(string(title_text), "_", "\_"), " ", "\ ");
title_handle = graph_thick_canvas_text(canvas_ax, ...
    "$\mathrm{" + title_text + "}$", title_font_size);
subtitle_handle = graph_thick_canvas_text(canvas_ax, ...
    graph_thick_option(options, 'subtitle', 'solid: without; dashed: with'), ...
    subtitle_font_size);
xlabel_handle = graph_thick_canvas_text(canvas_ax, ...
    'number of matrix-vector products', font_size);
note_handle = graph_thick_canvas_text(canvas_ax, ...
    'runtime (s): without / with', 80);
continuation_handle = graph_thick_canvas_text(canvas_ax, ...
    'adaptive solid curves continue beyond the shown range', 80);
canvas_texts = [title_handle, subtitle_handle, xlabel_handle, note_handle, continuation_handle];
drawnow;
ax.Position = axes_position;
drawnow;
inset = ax.TightInset;
text_extents = zeros(numel(canvas_texts), 4);
for text_index = 1:numel(canvas_texts)
    text_extents(text_index, :) = canvas_texts(text_index).Extent;
end
legend_size = lgd.Position(3:4);
margin = 80;
left_margin = max(500, inset(1) + margin);
right_margin = max(margin, inset(3) + margin);
figure_size(1) = max([figure_size(1), left_margin + panel_size(1) + right_margin, ...
    legend_size(1) + 2 * margin, max(text_extents(:, 3)) + 2 * margin]);
% Let the data panel use the width required by the typography and legend.
axes_position(1) = left_margin;
axes_position(3) = figure_size(1) - left_margin - right_margin;
continuation_bottom = margin;
note_bottom = continuation_bottom + text_extents(5, 4) + 20;
legend_bottom = note_bottom + text_extents(4, 4) + 35;
xlabel_bottom = legend_bottom + legend_size(2) + 45;
axes_position(2) = xlabel_bottom + text_extents(3, 4) + inset(2) + 50;
subtitle_bottom = axes_position(2) + panel_size(2) + inset(4) + 40;
title_bottom = subtitle_bottom + text_extents(2, 4) + 30;
figure_size(2) = max(figure_size(2), title_bottom + text_extents(1, 4) + margin);
fig.Position(3:4) = figure_size;
canvas_ax.Position = [0, 0, figure_size];
ax.Position = axes_position;
plot_center = axes_position(1) + axes_position(3) / 2;
title_handle.Position = [plot_center, title_bottom, 0];
subtitle_handle.Position = [plot_center, subtitle_bottom, 0];
xlabel_handle.Position = [plot_center, xlabel_bottom, 0];
note_handle.Position = [figure_size(1) / 2, note_bottom, 0];
continuation_handle.Position = [figure_size(1) / 2, continuation_bottom, 0];
lgd.Position = [(figure_size(1) - legend_size(1)) / 2, legend_bottom, legend_size];
drawnow;

% Include actual axes, text, and legend bounds in the exported canvas.
for layout_pass = 1:2
    inset = ax.TightInset;
    content_min = ax.Position(1:2) - inset(1:2);
    content_max = ax.Position(1:2) + ax.Position(3:4) + inset(3:4);
    for text_handle = canvas_texts
        extent = text_handle.Extent;
        content_min = min(content_min, extent(1:2));
        content_max = max(content_max, extent(1:2) + extent(3:4));
    end
    content_min = min(content_min, lgd.Position(1:2));
    content_max = max(content_max, lgd.Position(1:2) + lgd.Position(3:4));
    shift = max([0, 0], margin - content_min);
    position = ax.Position;
    position(1:2) = position(1:2) + shift;
    ax.Position = position;
    for text_handle = canvas_texts
        position = text_handle.Position;
        position(1:2) = position(1:2) + shift;
        text_handle.Position = position;
    end
    legend_position = lgd.Position;
    legend_position(1:2) = legend_position(1:2) + shift;
    lgd.Position = legend_position;
    figure_size = max(fig.Position(3:4), content_max + shift + margin);
    fig.Position(3:4) = figure_size;
    canvas_ax.Position = [0, 0, figure_size];
    drawnow;
end

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

function handle = graph_thick_canvas_text(ax, value, font_size)
handle = text(ax, 0, 0, value, 'Units', 'pixels', 'Position', [0, 0, 0], ...
    'Interpreter', 'latex', 'FontSize', font_size, 'FontWeight', 'normal', ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
    'Clipping', 'off');
end

function [limits, ticks] = graph_thick_zoom_range(maximum)
target_step = max(1, 1.1 * maximum / 4);
scale = 10 ^ floor(log10(target_step));
candidates = [1, 2, 2.5, 5, 10] * scale;
step = candidates(find(candidates >= target_step, 1));
limits = [0, 4 * step];
ticks = 0:step:limits(2);
end

function value = graph_thick_option(options, name, default_value)
if isfield(options, name)
    value = options.(name);
else
    value = default_value;
end
end
