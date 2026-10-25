exp_qcd_invsqrt_settings;

% Compare the average Arnoldi call cost for every method in each family.
selection.data_path = data_path;
selection.file_prefix = "qcd_invsqrt";
selection.krylov_dim = 150;
selection.thick = false;
group_names = ["Arnoldi", "HmArnoldi"];
method_groups = ["arnoldi", "sarnoldi", "last_orth"; ...
                 "hmarnoldi", "shmarnoldi", "last_hmorth"];
method_labels = ["Arn", "sArn", "adaArn"; ...
                 "HmArn", "sHmArn", "adaHmArn"];

for sfactor = sfactor_list
    selection.sfactor = sfactor;
    for group_index = 1:numel(group_names)
        baseline = qcd_read_median_run(selection, "fixed", ...
            method_groups(group_index, 1), inf);
        baseline_time = baseline.arnoldi_info.t_total / baseline.iter;
        result_list = struct('method_name', {}, 'group_name', {}, 'info', {});
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
                info = qcd_read_median_run(selection, mode, method, truncation_length);
                method_name = label;
                if mode == "adaptive"
                    method_name = label + "-" + string(truncation_length);
                end
                result_list(end + 1) = struct('method_name', method_name, ...
                    'group_name', group_names(group_index), 'info', info);
            end
        end
        if isempty(result_list)
            continue
        end

        figure_prefix = fullfile(figure_path, "qcd_invsqrt_arnoldi_ratios_" ...
            + lower(group_names(group_index)) ...
            + "_krylov_dim_" + string(selection.krylov_dim) ...
            + "_sfactor_" + string(sfactor));
        figure_title = "QCD";
        qcd_plot_ratios_arnoldi(result_list, figure_prefix, figure_title, ...
            baseline_time, method_labels(group_index, 1));
    end
end

% Local helpers keep QCD plotting independent of the other examples.
function info = qcd_read_median_run(selection, mode, method, truncation_length)
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

function qcd_plot_ratios_arnoldi(result_list, figure_prefix, figure_title, baseline_time, baseline_name)
% Compare fixed and adaptive methods using average time per Arnoldi call.
% Baselines are shared or supplied separately for each method.

if nargin < 5
    baseline_name = "Arn";
end

if isscalar(baseline_time)
    baseline_time = repmat(baseline_time, numel(result_list), 1);
else
    assert(numel(baseline_time) == numel(result_list), ...
        'Provide one baseline time per result, or a scalar baseline.');
    baseline_time = baseline_time(:);
end

if isempty(result_list)
    return
end

time_names = ["Apply operator", "Apply sketch", "Coefficients", ...
    "Subtraction", "Rank-1 update", "Other"];
ratios = zeros(numel(result_list), numel(time_names));
bar_labels = strings(1, numel(result_list));
group_names = strings(1, numel(result_list));
for it_result = 1 : numel(result_list)
    info = result_list(it_result).info;
    arnoldi_info = info.arnoldi_info;
    % The solver accumulates these times over info.iter Arnoldi calls.
    ratios(it_result, :) = [arnoldi_info.t_apply_op, arnoldi_info.t_apply_sketch, ...
        arnoldi_info.t_coeff, arnoldi_info.t_substract, ...
        arnoldi_info.t_update, arnoldi_info.t_other] / info.iter / baseline_time(it_result);
    bar_labels(it_result) = result_list(it_result).method_name;
    if isfield(result_list, 'group_name')
        group_names(it_result) = result_list(it_result).group_name;
    end
end

y_label = "run time relative to " + string(baseline_name);
qcd_plot_ratio_bars(ratios, time_names, bar_labels, figure_prefix, ...
    figure_title, y_label, group_names);

end

function qcd_plot_ratio_bars(ratios, time_names, bar_labels, figure_prefix, figure_title, y_label, group_names)
%PLOT_RATIO_BARS Draw timing parts relative to shared or per-group baselines.
% The input already contains the desired ratios; do not normalize each bar.
% Optional group_names inserts space and a divider between consecutive groups.

if isempty(bar_labels)
    return
end
if nargin < 7
    group_names = strings(1, numel(bar_labels));
else
    assert(numel(group_names) == numel(bar_labels), ...
        'Provide one group name per method.');
    group_names = reshape(string(group_names), 1, []);
end
has_groups = numel(unique(group_names(strlength(group_names) > 0))) > 1;
x_values = 1:numel(bar_labels);
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
active = any(ratios ~= 0, 1);
ratios = ratios(:, active);
time_names = time_names(active);
colors = colors(active, :);

font_size = 128;
legend_font_size = 104;
figure_size = [2300, 1650];
axes_position = [480, 400, 1100, 900];
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
    bars = bar(ax, x_values, ratios(:, stack_order), "stacked", ...
        "EdgeColor", "black", "LineWidth", 0.25);
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
maximum_height = max([1; sum(ratios, 2)]);
if has_groups
    ylim(ax, [0, 1.20 * maximum_height]);
else
    ylim(ax, [0, 1.05 * maximum_height]);
end
y_ticks = 0:0.2:ax.YLim(2);
yticks(ax, y_ticks);
yticklabels(ax, string(round(100 * y_ticks)) + "\%");
xticks(ax, x_values);
xticklabels(ax, replace(bar_labels, "_", "\_"));
xtickangle(ax, 40);
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
    for label_handle = [ax.XLabel, ax.YLabel, ax.Title]
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
