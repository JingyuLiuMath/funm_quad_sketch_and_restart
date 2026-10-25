function plot_ratio_bars(ratios, time_names, bar_labels, figure_prefix, figure_title, y_label, group_names, options)
%PLOT_RATIO_BARS Draw timing parts relative to shared or per-group baselines.
% The input already contains the desired ratios; do not normalize each bar.
% Optional group_names inserts space and a divider between consecutive groups.
% options.paired_bars places consecutive same-label results on either side
% of one method tick. options.subtitle describes their left/right order.

if isempty(bar_labels)
    return
end
if nargin < 8
    options = struct();
end
paired_bars = isfield(options, 'paired_bars') && options.paired_bars;
bar_labels = reshape(string(bar_labels), 1, []);
if nargin < 7
    group_names = strings(1, numel(bar_labels));
else
    assert(numel(group_names) == numel(bar_labels), ...
        'Provide one group name per method.');
    group_names = reshape(string(group_names), 1, []);
end
has_groups = numel(unique(group_names(strlength(group_names) > 0))) > 1;
if paired_bars
    assert(mod(numel(bar_labels), 2) == 0 && ...
        all(bar_labels(1:2:end) == bar_labels(2:2:end)), ...
        'Paired bars require consecutive results with matching method names.');
    assert(all(group_names(1:2:end) == group_names(2:2:end)), ...
        'The two bars in a method pair must belong to the same group.');
    method_centers = 1:(numel(bar_labels) / 2);
    x_values = reshape([method_centers - 0.18; method_centers + 0.18], 1, []);
else
    x_values = 1:numel(bar_labels);
end
group_starts = 1;
if has_groups
    boundaries = find(group_names(1:end - 1) ~= group_names(2:end));
    group_starts = [1, boundaries + 1];
    for boundary = boundaries
        x_values(boundary + 1:end) = x_values(boundary + 1:end) + 1.5;
    end
end

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
missing_bars = all(isnan(ratios), 2);
active = any(isfinite(ratios) & ratios ~= 0, 1);
ratios = ratios(:, active);
time_names = time_names(active);
colors = colors(active, :);

font_size = 128;
legend_font_size = 104;
figure_size = [2300, 1650];
axes_position = [480, 400, 1100, 900];
if paired_bars
    axes_position(3) = max([axes_position(3), 1900, ...
        440 * (numel(bar_labels) / 2)]);
    figure_size(1) = max(figure_size(1), axes_position(1) + axes_position(3) + 550);
end
if has_groups
    axes_position(3) = max(axes_position(3), 220 * numel(bar_labels));
    figure_size(1) = max(figure_size(1), axes_position(1) + axes_position(3) + 550);
end
fig = figure("Color", "white", "Units", "pixels", ...
    "Position", [100, 100, figure_size], "Renderer", "painters");
ax = axes(fig, "Units", "pixels", "Position", axes_position);
lgd = gobjects(0);

% Stack the first timing part at the top while retaining the legend order.
stack_order = numel(time_names):-1:1;
if ~isempty(stack_order)
    % MATLAB scales stacked bar width by the smallest spacing in XData.
    % For paired positions +/-0.18, width 0.8 leaves a small within-pair gap.
    bar_width = 0.8;
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
if has_groups
    ylim(ax, [0, 1.20 * maximum_height]);
else
    ylim(ax, [0, 1.05 * maximum_height]);
end
if paired_bars
    target_step = ax.YLim(2) / 4;
    tick_scale = 10 ^ floor(log10(target_step));
    tick_steps = [1, 2, 2.5, 5, 10] * tick_scale;
    tick_step = tick_steps(find(tick_steps >= target_step, 1));
    upper_tick = ceil(ax.YLim(2) / tick_step) * tick_step;
    ylim(ax, [0, upper_tick]);
    y_ticks = 0:tick_step:upper_tick;
else
    y_ticks = 0:0.2:ax.YLim(2);
end
yticks(ax, y_ticks);
yticklabels(ax, string(round(100 * y_ticks)) + "\%");
if paired_bars
    xticks(ax, mean(reshape(x_values, 2, []), 1));
    xticklabels(ax, replace(bar_labels(1:2:end), "_", "\_"));
else
    xticks(ax, x_values);
    xticklabels(ax, replace(bar_labels, "_", "\_"));
end
xtickangle(ax, 40);
if paired_bars && any(missing_bars)
    text(ax, x_values(missing_bars), ...
        repmat(0.025 * ax.YLim(2), 1, nnz(missing_bars)), "N/A", ...
        "Interpreter", "latex", "FontSize", legend_font_size, ...
        "FontWeight", "normal", "HorizontalAlignment", "center", ...
        "VerticalAlignment", "bottom", "HandleVisibility", "off");
end
if has_groups
    group_ends = [group_starts(2:end) - 1, numel(bar_labels)];
    hold(ax, "on");
    for group_index = 1:numel(group_starts)
        group_center = mean(x_values([group_starts(group_index), group_ends(group_index)]));
        text(ax, group_center, 0.98 * ax.YLim(2), ...
            replace(group_names(group_starts(group_index)), "_", "\_"), ...
            "Interpreter", "latex", "FontSize", font_size, ...
            "HorizontalAlignment", "center", "VerticalAlignment", "top", ...
            "HandleVisibility", "off");
        if group_index < numel(group_starts)
            divider = mean(x_values([group_ends(group_index), group_starts(group_index + 1)]));
            xline(ax, divider, "--", "Color", [0.4, 0.4, 0.4], ...
                "LineWidth", 1, "HandleVisibility", "off");
        end
    end
    hold(ax, "off");
end
ylabel(ax, y_label, "Interpreter", "latex", "FontSize", font_size);
title_text = replace(replace(string(figure_title), "_", "\_"), " ", "\ ");
title(ax, "$\mathrm{" + title_text + "}$", "Interpreter", "latex", ...
    "FontSize", font_size, "FontWeight", "normal");
subtitle_handle = gobjects(0);
if isfield(options, 'subtitle') && strlength(string(options.subtitle)) > 0
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
