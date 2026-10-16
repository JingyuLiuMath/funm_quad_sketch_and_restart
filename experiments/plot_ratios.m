function plot_ratios(result_list, figure_prefix, figure_title)
% Each entry contains the label and info from one selected run.

time_names = ["Arnoldi", "Initial action", "Quadrature", ...
    "Thick restart", "Update", "Other"];
fractions = zeros(numel(result_list), numel(time_names));
bar_labels = strings(1, numel(result_list));
for it_result = 1 : numel(result_list)
    info = result_list(it_result).info;
    fractions(it_result, :) = [info.t_arnoldi, info.t_funm_init, ...
        info.t_funm_quad, info.t_thick, info.t_update, info.t_other] / info.t_total;
    bar_labels(it_result) = result_list(it_result).method_name;
end

fig = figure("Color", "white", "Position", [100, 100, 1200, 750]);
ax = axes(fig);
bar(ax, fractions, "stacked");
grid(ax, "on");
ylim(ax, [0, 1]);
xticks(ax, 1 : numel(bar_labels));
xticklabels(ax, bar_labels);
xtickangle(ax, 30);
set(ax, "TickLabelInterpreter", "none");
ylabel(ax, "Fraction");
title(ax, figure_title, "Interpreter", "none");
legend(ax, time_names, "Interpreter", "none", "Orientation", "horizontal", ...
    "NumColumns", 3, "Location", "southoutside");
exportgraphics(fig, figure_prefix + ".pdf", "ContentType", "vector");
exportgraphics(fig, figure_prefix + ".png", "Resolution", 300);

end
