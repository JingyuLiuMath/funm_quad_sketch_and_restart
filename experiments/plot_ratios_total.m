function plot_ratios_total(result_list, figure_prefix, figure_title, baseline_time)
% Normalize timing parts by the ordinary fixed Arnoldi total time.

time_names = ["Arnoldi", "Initial action", "Quadrature", ...
    "Thick restart", "Update", "Other"];
ratios = zeros(numel(result_list), numel(time_names));
bar_labels = strings(1, numel(result_list));
for it_result = 1 : numel(result_list)
    info = result_list(it_result).info;
    ratios(it_result, :) = [info.t_arnoldi, info.t_funm_init, ...
        info.t_funm_quad, info.t_thick, info.t_update, info.t_other] / baseline_time;
    bar_labels(it_result) = result_list(it_result).method_name;
end

fig = figure("Color", "white", "Position", [100, 100, 1200, 750]);
ax = axes(fig);
bar(ax, ratios, "stacked");
grid(ax, "on");
ylim(ax, [0, 1.05 * max([1; sum(ratios, 2)])]);
y_ticks = 0 : 0.2 : ax.YLim(2);
yticks(ax, y_ticks);
yticklabels(ax, string(round(100 * y_ticks)) + "%");
xticks(ax, 1 : numel(bar_labels));
xticklabels(ax, bar_labels);
xtickangle(ax, 30);
set(ax, "TickLabelInterpreter", "none");
ylabel(ax, "Run time relative to Arnoldi");
title(ax, figure_title, "Interpreter", "none");
legend(ax, time_names, "Interpreter", "none", "Orientation", "horizontal", ...
    "NumColumns", 3, "Location", "southoutside");
exportgraphics(fig, figure_prefix + ".pdf", "ContentType", "vector");
exportgraphics(fig, figure_prefix + ".png", "Resolution", 300);

end
