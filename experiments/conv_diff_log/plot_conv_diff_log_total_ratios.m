exp_conv_diff_log_settings;

% Fixed methods are paired under one tick, ordered left m=150, right m=200.
fixed_krylov_dims = [150, 200];
assert(all(ismember(fixed_krylov_dims, krylov_dim_list)), ...
    'The ratio plots require fixed results at m=150 and m=200.');
selection.data_path = data_path;
selection.file_prefix = "conv_diff_log";
selection.thick = false;
group_names = ["Arnoldi", "HmArnoldi"];
method_groups = ["arnoldi", "sarnoldi", "last_orth"; ...
                 "hmarnoldi", "shmarnoldi", "last_hmorth"];
method_labels = ["Arn", "sArn", "adaArn"; ...
                 "HmArn", "sHmArn", "adaHmArn"];
plot_options.subtitle = ["Fixed $m$: $150$ (left), $200$ (right)"; ...
    "Adaptive: $m_{\max}=" + string(adaptive_krylov_dim) + "$"];

for sfactor = sfactor_list
    selection.sfactor = sfactor;
    % Share the same Arn(m=200) reference across both dimensions and families.
    selection.krylov_dim = 200;
    baseline = conv_read_median_run(selection, "fixed", "arnoldi", inf);
    baseline_time = baseline.t_total;
    baseline_name = "Arn";
    for group_index = 1:numel(group_names)
        result_list = struct('method_name', {}, 'info', {});
        for method_index = 1:2
            method = method_groups(group_index, method_index);
            if ~ismember(method, fixed_method_list)
                continue
            end
            for krylov_dim = fixed_krylov_dims
                selection.krylov_dim = krylov_dim;
                info = conv_read_median_run(selection, "fixed", method, inf);
                result_list(end + 1) = struct('method_name', method_labels(group_index, method_index), ...
                    'info', info); %#ok<SAGROW>
            end
        end
        method = method_groups(group_index, 3);
        if ismember(method, adaptive_update_list)
            selection.krylov_dim = adaptive_krylov_dim;
            for truncation_length = truncation_length_list
                info = conv_read_median_run(selection, "adaptive", method, truncation_length);
                label = method_labels(group_index, 3) + "-" + string(truncation_length);
                result_list(end + 1) = struct('method_name', label, ...
                    'info', info); %#ok<SAGROW>
            end
        end
        if isempty(result_list)
            continue
        end

        figure_prefix = fullfile(figure_path, "conv_diff_log_total_ratios_" ...
            + lower(group_names(group_index)) + "_krylov_dim_200" ...
            + "_sfactor_" + string(sfactor));
        conv_plot_ratios_total(result_list, figure_prefix, "Conv-Diff", ...
            baseline_time, baseline_name, plot_options);
    end
end

% Local helpers keep this convection-diffusion plot independent of other examples.
function info = conv_read_median_run(selection, mode, method, truncation_length)
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

function conv_plot_ratios_total(result_list, figure_prefix, figure_title, baseline_time, baseline_name, options)
% Normalize timing parts by a shared or per-method baseline total time.
% Optional baseline_name labels the reference method; options is forwarded
% to conv_plot_ratio_bars (including subtitle and bar_width).

if nargin < 5
    baseline_name = "Arn";
end
if nargin < 6
    options = struct();
end

if isscalar(baseline_time)
    baseline_time = repmat(baseline_time, numel(result_list), 1);
else
    assert(numel(baseline_time) == numel(result_list), ...
        'Provide one baseline time per result, or a scalar baseline.');
    baseline_time = baseline_time(:);
end

time_names = ["Basis generation", "Initial action", "Quadrature", ...
    "Thick restart", "Update", "Other"];
ratios = zeros(numel(result_list), numel(time_names));
bar_labels = strings(1, numel(result_list));
for it_result = 1 : numel(result_list)
    info = result_list(it_result).info;
    if isempty(info)
        % An unsupported method variant has no measured timing breakdown.
        ratios(it_result, :) = NaN;
    else
        ratios(it_result, :) = [info.t_arnoldi, info.t_funm_init, ...
            info.t_funm_quad, info.t_thick, info.t_update, info.t_other] / baseline_time(it_result);
    end
    bar_labels(it_result) = result_list(it_result).method_name;
end

y_label = "run time / $t_{\mathrm{" + string(baseline_name) + "},200}$";
conv_plot_ratio_bars(ratios, time_names, bar_labels, figure_prefix, ...
    figure_title, y_label, options);

end

function conv_plot_ratio_bars(ratios, time_names, bar_labels, figure_prefix, figure_title, y_label, options)
%PLOT_RATIO_BARS Draw stacked timing parts with one tick per method.
% Consecutive fixed results share a label: m=150 is left, m=200 is right.
% Adaptive labels occur once and receive a single centered bar.

if isempty(bar_labels)
    return
end
if nargin < 7
    options = struct();
end
bar_labels = reshape(string(bar_labels), 1, []);
assert(size(ratios, 1) == numel(bar_labels), ...
    'Provide one timing row per bar.');
label_starts = [1, find(bar_labels(1:end - 1) ~= bar_labels(2:end)) + 1];
label_ends = [label_starts(2:end) - 1, numel(bar_labels)];
assert(all(label_ends - label_starts + 1 <= 2), ...
    'Each method has one bar or a fixed-dimension pair.');
% Fixed pairs touch; distinct methods have a quarter-bar-width gap,
% matching the spacing-to-width ratio of the other examples' BarWidth=0.8.
method_boundaries = [0, bar_labels(1:end - 1) ~= bar_labels(2:end)];
x_values = (1:numel(bar_labels)) + 0.25 * cumsum(method_boundaries);
method_centers = (x_values(label_starts) + x_values(label_ends)) / 2;

palette = [0.35, 0.66, 0.82; ...
           0.89, 0.43, 0.37; ...
           0.38, 0.70, 0.47; ...
           0.95, 0.69, 0.29; ...
           0.44, 0.57, 0.75; ...
           0.68, 0.52, 0.72];
% Keep timing meanings tied to their colors when a component is omitted.
color_names = ["Arnoldi", "Basis generation", "Initial action", "Quadrature", "Thick restart", ...
    "Update", "Rank-1 update", "Other", "Apply operator", "Apply sketch", ...
    "Coefficients", "Subtraction"];
color_indices = [1, 1, 2, 3, 4, 5, 5, 6, 1, 2, 3, 4];
time_names = reshape(string(time_names), 1, []);
colors = palette(1 + mod(0:numel(time_names) - 1, size(palette, 1)), :);
[known_names, name_indices] = ismember(time_names, color_names);
colors(known_names, :) = palette(color_indices(name_indices(known_names)), :);
active = any(isfinite(ratios) & ratios ~= 0, 1);
ratios = ratios(:, active);
time_names = time_names(active);
colors = colors(active, :);

font_size = 128;
legend_font_size = 104;
% Match the reference plotting area and PDF width for equal-width LaTeX inclusion.
axes_position = [480, 400, 1100, 900];
figure_size = [2544, 1650];
fig = figure("Color", "white", "Units", "pixels", ...
    "Position", [100, 100, figure_size], "Renderer", "painters");
ax = axes(fig, "Units", "pixels", "Position", axes_position);
lgd = gobjects(0);

% Stack the first timing part at the top while retaining the legend order.
stack_order = numel(time_names):-1:1;
if ~isempty(stack_order)
    % Keep an actual width of one even if only adaptive single bars remain.
    bar_width = 1;
    if numel(x_values) > 1
        bar_width = 1 / min(diff(x_values));
    end
    if isfield(options, 'bar_width')
        bar_width = options.bar_width;
    end
    bars = bar(ax, x_values, ratios(:, stack_order), "stacked", ...
        "BarWidth", bar_width, "EdgeColor", "black", "LineWidth", 0.25);
    for index = 1:numel(bars)
        bars(index).FaceColor = colors(stack_order(index), :);
    end
    lgd = legend(ax, bars(end:-1:1), replace(time_names, "_", "\_"), ...
        "Interpreter", "latex", "FontSize", legend_font_size, ...
        "Location", "eastoutside", "NumColumns", 1, "Box", "off");
    lgd.ItemTokenSize = [80, 60];
end

set(ax, "FontSize", font_size, "FontWeight", "normal", ...
    "PositionConstraint", "innerposition", "TickLabelInterpreter", "latex", ...
    "XGrid", "off", "YGrid", "on", "GridLineStyle", "--", ...
    "Box", "off", "Layer", "top");
bar_heights = sum(ratios, 2);
maximum_height = max([1; bar_heights(isfinite(bar_heights))]);
ylim(ax, [0, 1.05 * maximum_height]);
% Limit tick density for every bar layout, including ratios well above 100%.
target_step = ax.YLim(2) / 5;
tick_scale = 10 ^ floor(log10(target_step));
tick_steps = [1, 2, 2.5, 5, 10] * tick_scale;
tick_step = tick_steps(find(tick_steps >= target_step, 1));
upper_tick = ceil(ax.YLim(2) / tick_step) * tick_step;
ylim(ax, [0, upper_tick]);
y_ticks = 0:tick_step:upper_tick;
yticks(ax, y_ticks);
yticklabels(ax, string(round(100 * y_ticks)) + "\%");
xticks(ax, method_centers);
xticklabels(ax, replace(bar_labels(label_starts), "_", "\_"));
% Leave room between the three adaptive labels at the reference font size.
xtickangle(ax, 50);
xlim(ax, [x_values(1) - 0.75, x_values(end) + 0.75]);
ylabel(ax, y_label, "Interpreter", "latex", "FontSize", font_size);
title_text = replace(replace(string(figure_title), "_", "\_"), " ", "\ ");
title(ax, "$\mathrm{" + title_text + "}$", "Interpreter", "latex", ...
    "FontSize", font_size, "FontWeight", "normal");
subtitle_handle = gobjects(0);
if isfield(options, 'subtitle') && any(strlength(string(options.subtitle)) > 0)
    subtitle_handle = subtitle(ax, options.subtitle, "Interpreter", "latex", ...
        "FontSize", legend_font_size, "FontWeight", "normal");
end
drawnow;

% Preserve the plot box and leave room for long labels at the shared font size.
set(ax, "Position", axes_position);
drawnow;
insets = ax.TightInset;
axes_position(1:2) = axes_position(1:2) + ...
    max([0, 0], insets(1:2) + 25 - axes_position(1:2));
ax.Position = axes_position;
figure_size = max(figure_size, axes_position(1:2) + axes_position(3:4) + ...
    insets(3:4) + 25);
if ~isempty(lgd)
    set(lgd, "Units", "pixels");
    legend_position = lgd.Position;
    legend_position(1) = axes_position(1) + axes_position(3) + 25;
    legend_position(2) = axes_position(2) + (axes_position(4) - legend_position(4)) / 2;
    lgd.Position = legend_position;
    figure_size(1) = max(figure_size(1), legend_position(1) + legend_position(3) + 40);
end
figure_position = fig.Position;
figure_position(3:4) = figure_size;
fig.Position = figure_position;
drawnow;

% TightInset can miss the upper end of a long rotated ylabel. Check each
% label's actual pixel rectangle after layout, without enlarging the plot box.
safety_margin = 50;
for layout_pass = 1:2
    axes_position = ax.Position;
    insets = ax.TightInset;
    content_min = axes_position(1:2) - insets(1:2);
    content_max = axes_position(1:2) + axes_position(3:4) + insets(3:4);
    label_handles = [ax.XLabel, ax.YLabel, ax.Title];
    if ~isempty(subtitle_handle)
        label_handles = [label_handles, subtitle_handle];
    end
    for label_handle = label_handles
        original_units = label_handle.Units;
        label_handle.Units = "pixels";
        label_extent = label_handle.Extent;
        label_handle.Units = original_units;
        label_min = axes_position(1:2) + label_extent(1:2);
        label_max = label_min + label_extent(3:4);
        content_min = min(content_min, label_min);
        content_max = max(content_max, label_max);
    end
    if ~isempty(lgd)
        legend_position = lgd.Position;
        content_min = min(content_min, legend_position(1:2));
        content_max = max(content_max, legend_position(1:2) + legend_position(3:4));
    end

    layout_shift = max([0, 0], safety_margin - content_min);
    axes_position(1:2) = axes_position(1:2) + layout_shift;
    ax.Position = axes_position;
    if ~isempty(lgd)
        legend_position(1:2) = legend_position(1:2) + layout_shift;
        lgd.Position = legend_position;
    end
    figure_position = fig.Position;
    figure_position(3:4) = max(figure_position(3:4), ...
        content_max + layout_shift + safety_margin);
    fig.Position = figure_position;
    drawnow;
end

figure_prefix = string(figure_prefix);
figure_directory = fileparts(figure_prefix);
if strlength(figure_directory) > 0 && ~isfolder(figure_directory)
    mkdir(figure_directory);
end
% Use the same full-canvas export as the error-history figures on R2023b.
paper_size = fig.Position(3:4) / get(groot, "ScreenPixelsPerInch");
set(fig, "PaperUnits", "inches", "PaperSize", paper_size, ...
    "PaperPosition", [0, 0, paper_size], "PaperPositionMode", "manual");
print(fig, char(figure_prefix + ".pdf"), "-dpdf", "-painters");
print(fig, char(figure_prefix + ".png"), "-dpng", "-r300");

end
