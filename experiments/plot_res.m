function plot_res(result_list, figure_prefix, figure_title)
% Each entry contains a selected history and its method/truncation indices.

fig = figure("Color", "white", "Position", [100, 100, 1100, 800]);
ax = axes(fig);
method_colors = lines(max([result_list.method_index]));
line_styles = ["-", "--", ":", "-."];
markers = ["o", "x", "s", "d"];
for it_result = 1 : numel(result_list)
    result = result_list(it_result);
    it_style = 1 + mod(result.truncation_index - 1, numel(line_styles));
    semilogy(ax, 1 : numel(result.res_hist), result.res_hist, ...
        "Color", method_colors(result.method_index, :), ...
        "LineStyle", line_styles(it_style), "Marker", markers(it_style), ...
        "LineWidth", 2, "MarkerSize", 4, "DisplayName", result.method_name);
    hold(ax, "on");
end
grid(ax, "on");
xlabel(ax, "Restart");
ylabel(ax, "True relative error");
title(ax, figure_title, "Interpreter", "none");
legend(ax, "Interpreter", "none", "Location", "best", "NumColumns", 2);
exportgraphics(fig, figure_prefix + ".pdf", "ContentType", "vector");
exportgraphics(fig, figure_prefix + ".png", "Resolution", 300);

end
