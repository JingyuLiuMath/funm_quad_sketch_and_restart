exp_graph_laplacian_invsqrt_settings;

% Select the same runs and methods as plot_graph_laplacian_invsqrt_err.
selection.data_path = data_path;
selection.file_prefix = "frac_laplacian_invsqrt";
selection.krylov_dim = 150;
selection.thick = false;
group_names = ["arnoldi", "hmarnoldi"];
method_groups = ["arnoldi", "sarnoldi", "last_orth"; ...
                 "hmarnoldi", "shmarnoldi", "last_hmorth"];
method_labels = ["Arn", "sArn", "adaArn"; ...
                 "HmArn", "sHmArn", "adaHmArn"];

plot_options.figure_size = [2300, 1650];
plot_options.axes_position = [480, 400, 1100, 900];
plot_options.font_size = 128;
plot_options.legend_font_size = 104;
plot_options.line_width = 12;
plot_options.marker_size = 40;
plot_options.legend_location = "eastoutside";
plot_options.legend_gap = 25;
plot_options.legend_token_size = [80, 60];
plot_options.x_label = "iteration";
plot_options.y_label = "quadrature points";

for sfactor = sfactor_list
    selection.sfactor = sfactor;
    for group_index = 1:numel(group_names)
        result_list = struct('method_name', {}, 'num_quadpoints', {});
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
                info = graph_read_median_run(selection, mode, method, truncation_length);
                assert(numel(info.num_quadpoints) == info.iter, ...
                    'The quadrature history must contain one entry per iteration.');
                result.method_name = label;
                if mode == "adaptive"
                    result.method_name = label + "-" + string(truncation_length);
                end
                % The first cycle uses initial_action, so its count is zero.
                % Later counts are the final fine-rule sizes, not cumulative work.
                result.num_quadpoints = info.num_quadpoints(:).';
                result_list(end + 1) = result;
            end
        end
        if isempty(result_list)
            continue
        end

        % Use four readable iteration ticks, as in the error plots.
        max_iter = max(arrayfun(@(result) numel(result.num_quadpoints), result_list));
        target_step = max(1, 1.03 * max_iter / 3);
        step_scale = 10 ^ floor(log10(target_step));
        tick_step = ceil(target_step / step_scale) * step_scale;
        plot_options.x_limits = [0, 3 * tick_step];
        plot_options.x_ticks = 0:tick_step:3 * tick_step;

        % Quadrature refinement doubles the rule size; keep zero on a linear axis.
        max_points = max([result_list.num_quadpoints]);
        quad_tick_step = 2 ^ max(0, ceil(log2(max(1, max_points) / 4)));
        quad_upper_tick = max(quad_tick_step, ...
            ceil(max_points / quad_tick_step) * quad_tick_step);
        plot_options.y_limits = [0, 1.1 * quad_upper_tick];
        plot_options.y_ticks = 0:quad_tick_step:quad_upper_tick;

        figure_prefix = fullfile(figure_path, "graph_laplacian_invsqrt_quadpoints_" ...
            + group_names(group_index) + "_sfactor_" + string(sfactor));
        graph_plot_quadpoints(result_list, figure_prefix, "Graph Laplacian", plot_options);
    end
end

% Local helpers keep this Graph Laplacian plot independent of other examples.
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

function fig = graph_plot_quadpoints(result_list, figure_prefix, matname, options)
%PLOT_QUADPOINTS Plot quadrature-point counts for every restart cycle.
% Each result contains method_name and num_quadpoints. The full history,
% including the initial zero count, is plotted against 1:numel(history).
% options supports the presentation settings in graph_plot_res, with x/y limits
% and ticks. Export a PDF and PNG using figure_prefix.

if nargin < 4
    options = struct;
end
figure_size = graph_plot_option(options, "figure_size", [1200, 850]);
legend_location = graph_plot_option(options, "legend_location", "northeast");
legend_gap = graph_plot_option(options, "legend_gap", 50);
legend_right_padding = graph_plot_option(options, "legend_right_padding", 40);
x_label = graph_plot_option(options, "x_label", "iteration");
y_label = graph_plot_option(options, "y_label", "number of quadrature points");
axes_position = graph_plot_option(options, "axes_position", [250, 190, 847, 563]);
font_size = graph_plot_option(options, "font_size", 80);
legend_font_size = graph_plot_option(options, "legend_font_size", 0.8 * font_size);
line_width = graph_plot_option(options, "line_width", 6);
marker_size = graph_plot_option(options, "marker_size", 20);
line_styles = ["-", "--", "-.", ":"];
markers = ["o", "s", "^", "d", "v", ">", "<", "p"];
method_colors = lines(numel(result_list));

fig = figure("Color", "white", "Units", "pixels", ...
    "Position", [100, 100, figure_size], "Renderer", "painters");
ax = axes(fig, "Units", "pixels", "Position", axes_position, ...
    "PositionConstraint", "innerposition", "FontSize", font_size, ...
    "FontWeight", "normal", "TickLabelInterpreter", "latex", "Box", "on");
hold(ax, "on");
all_counts = [];

for it_method = 1:numel(result_list)
    result = result_list(it_method);
    assert(all(isfield(result, {'method_name', 'num_quadpoints'})), ...
        'Each result must contain method_name and num_quadpoints.');
    counts = result.num_quadpoints;
    assert(isnumeric(counts) && isreal(counts) && ...
        (isempty(counts) || isvector(counts)) && ...
        all(isfinite(counts(:))) && all(counts(:) >= 0) && ...
        all(counts(:) == fix(counts(:))), ...
        'num_quadpoints must be a vector of finite nonnegative integers.');
    counts = counts(:);
    iterations = (1:numel(counts)).';

    marker_step = max(1, ceil(numel(counts) / 12));
    marker_start = 1 + floor((it_method - 1) * marker_step / numel(result_list));
    line_style = line_styles(1 + mod(it_method - 1, numel(line_styles)));
    marker = markers(1 + mod(it_method - 1, numel(markers)));
    display_name = replace(string(result.method_name), "_", "\_");
    plot(ax, iterations, counts, "Color", method_colors(it_method, :), ...
        "LineStyle", line_style, "LineWidth", line_width, "Marker", marker, ...
        "MarkerSize", marker_size, "MarkerFaceColor", "none", ...
        "MarkerIndices", marker_start:marker_step:numel(counts), ...
        "DisplayName", display_name);
    all_counts = [all_counts; double(counts)]; %#ok<AGROW>
end

set(ax, "XScale", "linear", "YScale", "linear");
grid(ax, "off");
if isempty(all_counts)
    y_limits = [0, 1];
else
    y_limits = [0, max(1, 1.05 * max(all_counts))];
end
ylim(ax, graph_plot_option(options, "y_limits", y_limits));
if isfield(options, "x_limits")
    xlim(ax, options.x_limits);
end
if isfield(options, "x_ticks")
    xticks(ax, options.x_ticks);
end
if isfield(options, "y_ticks")
    yticks(ax, options.y_ticks);
end
xtickangle(ax, 0);
xlabel(ax, x_label, "Interpreter", "latex", "FontSize", font_size, ...
    "FontWeight", "normal");
ylabel(ax, y_label, "Interpreter", "latex", "FontSize", font_size, ...
    "FontWeight", "normal");
title_text = replace(replace(string(matname), "_", "\_"), " ", "\ ");
title(ax, "$\mathrm{" + title_text + "}$", ...
    "Interpreter", "latex", "FontSize", font_size, "FontWeight", "normal");
lgd = legend(ax, "Interpreter", "latex", "FontSize", legend_font_size, ...
    "FontWeight", "normal", "Location", legend_location);
if isfield(options, "legend_token_size")
    lgd.ItemTokenSize = options.legend_token_size;
end

% Add the margins needed by the chosen labels without reducing the axes.
drawnow;
set(ax, "Position", axes_position);
drawnow;
insets = ax.TightInset;
axes_position(1:2) = max(axes_position(1:2), insets(1:2) + 25);
ax.Position = axes_position;
figure_position = fig.Position;
figure_position(3:4) = max(figure_position(3:4), ...
    axes_position(1:2) + axes_position(3:4) + insets(3:4) + 25);
fig.Position = figure_position;

if strcmpi(legend_location, "eastoutside")
    drawnow;
    set(ax, "Units", "pixels", "Position", axes_position);
    set(lgd, "Units", "pixels");
    legend_position = lgd.Position;
    legend_position(1) = axes_position(1) + axes_position(3) + legend_gap;
    legend_position(2) = max(25, axes_position(2) + ...
        (axes_position(4) - legend_position(4)) / 2);
    lgd.Position = legend_position;
    figure_position = fig.Position;
    figure_position(3) = max(figure_position(3), ...
        legend_position(1) + legend_position(3) + legend_right_padding);
    figure_position(4) = max(figure_position(4), ...
        legend_position(2) + legend_position(4) + 25);
    fig.Position = figure_position;
    drawnow;
    ax.Position = axes_position;
    lgd.Position = legend_position;
end
drawnow;

% Export the full canvas, as in graph_plot_res on MATLAB R2023b.
figure_path = fileparts(figure_prefix);
if strlength(string(figure_path)) > 0 && ~isfolder(figure_path)
    mkdir(figure_path);
end
pixels_per_inch = get(groot, "ScreenPixelsPerInch");
paper_size = fig.Position(3:4) / pixels_per_inch;
set(fig, "PaperUnits", "inches", "PaperSize", paper_size, ...
    "PaperPosition", [0, 0, paper_size], "PaperPositionMode", "manual");
print(fig, char(string(figure_prefix) + ".pdf"), "-dpdf", "-painters");
print(fig, char(string(figure_prefix) + ".png"), "-dpng", "-r300");
end

function value = graph_plot_option(options, name, default_value)
if isfield(options, name)
    value = options.(name);
else
    value = default_value;
end
end
