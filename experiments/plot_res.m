function fig = plot_res(result_list, figure_prefix, matname, options)
%PLOT_RES Plot selected error histories and export a PDF and PNG.
% Each result contains method_name and res_hist. Optional matvec_hist gives
% cumulative matrix-vector counts; otherwise the abscissa is the iteration.
% options may set figure_size, axes_position, legend_location, legend_gap,
% legend_right_padding, legend_token_size, font_size, legend_font_size,
% line_width, marker_size, x_label, y_label, x_limits, x_ticks, and y_limits.
% Three-argument callers use the standard presentation settings.

if nargin < 4
    options = struct;
end
figure_size = option_value(options, "figure_size", [1200, 850]);
legend_location = option_value(options, "legend_location", "northeast");
legend_gap = option_value(options, "legend_gap", 50);
legend_right_padding = option_value(options, "legend_right_padding", 40);
x_label = option_value(options, "x_label", "iteration");
y_label = option_value(options, "y_label", "relative error");

axes_position = option_value(options, "axes_position", [250, 190, 847, 563]);
font_size = option_value(options, "font_size", 80);
legend_font_size = option_value(options, "legend_font_size", 0.8 * font_size);
line_width = option_value(options, "line_width", 6);
marker_size = option_value(options, "marker_size", 20);
line_styles = ["-", "--", "-.", ":"];
markers = ["o", "s", "^", "d", "v", ">", "<", "p"];
method_colors = lines(numel(result_list));

fig = figure("Color", "white", "Units", "pixels", ...
    "Position", [100, 100, figure_size], "Renderer", "painters");
ax = axes(fig, "Units", "pixels", "Position", axes_position, ...
    "PositionConstraint", "innerposition", "FontSize", font_size, ...
    "TickLabelInterpreter", "latex", "Box", "on");
hold(ax, "on");
positive_errors = [];

for it_method = 1:numel(result_list)
    result = result_list(it_method);
    errors = result.res_hist(:);
    if isfield(result, "matvec_hist") && ~isempty(result.matvec_hist)
        x_values = result.matvec_hist(:);
    else
        x_values = (1:numel(errors)).';
    end
    assert(numel(x_values) == numel(errors), ...
        'The error history and abscissa must have the same length.');

    marker_step = max(1, ceil(numel(errors) / 12));
    marker_start = 1 + floor((it_method - 1) * marker_step / numel(result_list));
    line_style = line_styles(1 + mod(it_method - 1, numel(line_styles)));
    marker = markers(1 + mod(it_method - 1, numel(markers)));
    display_name = replace(string(result.method_name), "_", "\_");
    semilogy(ax, x_values, errors, "Color", method_colors(it_method, :), ...
        "LineStyle", line_style, "LineWidth", line_width, "Marker", marker, ...
        "MarkerSize", marker_size, "MarkerFaceColor", "none", ...
        "MarkerIndices", marker_start:marker_step:numel(errors), ...
        "DisplayName", display_name);
    positive_errors = [positive_errors; errors(isfinite(errors) & errors > 0)]; %#ok<AGROW>
end

% Keep the reference decade ticks and extend them to include every error.
lower_exponent = -10;
upper_exponent = 0;
if ~isempty(positive_errors)
    lower_exponent = min(lower_exponent, 2 * floor(log10(min(positive_errors)) / 2));
    upper_exponent = max(upper_exponent, 2 * ceil(log10(max(positive_errors)) / 2));
end
set(ax, "YScale", "log", "YTick", 10 .^ (lower_exponent:2:upper_exponent), ...
    "YLim", [10 ^ lower_exponent, 10 ^ upper_exponent]);
grid(ax, "off");
if isfield(options, "y_limits")
    ylim(ax, options.y_limits);
end
if isfield(options, "x_limits")
    xlim(ax, options.x_limits);
end
if isfield(options, "x_ticks")
    xticks(ax, options.x_ticks);
end
xtickangle(ax, 0);
xlabel(ax, x_label, "Interpreter", "latex", "FontSize", font_size);
ylabel(ax, y_label, "Interpreter", "latex", "FontSize", font_size);
title_text = replace(replace(string(matname), "_", "\_"), " ", "\ ");
title(ax, "$\mathrm{" + title_text + "}$", ...
    "Interpreter", "latex", "FontSize", font_size, "FontWeight", "normal");
lgd = legend(ax, "Interpreter", "latex", "FontSize", legend_font_size, ...
    "Location", legend_location);
if isfield(options, "legend_token_size")
    lgd.ItemTokenSize = options.legend_token_size;
end

% Add only the margin needed by the chosen labels and font sizes.
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

% An outside legend must not reduce the requested plotting area.
if strcmpi(legend_location, "eastoutside")
    drawnow;
    set(ax, "Units", "pixels", "Position", axes_position);
    set(lgd, "Units", "pixels");
    legend_position = lgd.Position;
    legend_position(1) = axes_position(1) + axes_position(3) + legend_gap;
    legend_position(2) = axes_position(2) + (axes_position(4) - legend_position(4)) / 2;
    lgd.Position = legend_position;
    figure_position = fig.Position;
    figure_position(3) = max(figure_position(3), ...
        legend_position(1) + legend_position(3) + legend_right_padding);
    fig.Position = figure_position;
end
drawnow;

% print supports this full-canvas export on MATLAB R2023b.
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

function value = option_value(options, name, default_value)
if isfield(options, name)
    value = options.(name);
else
    value = default_value;
end
end
