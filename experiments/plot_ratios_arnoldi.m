function plot_ratios_arnoldi(result_list, figure_prefix, figure_title, baseline_time)
% Fixed cases only; baseline_time is the ordinary fixed Arnoldi time per call.

is_fixed = arrayfun(@(result) result.info.param.mode == "fixed", result_list);
result_list = result_list(is_fixed);
if isempty(result_list)
    return
end

time_names = ["Apply operator", "Apply sketch", "Coefficients", ...
    "Subtraction", "Update", "Other"];
ratios = zeros(numel(result_list), numel(time_names));
bar_labels = strings(1, numel(result_list));
for it_result = 1 : numel(result_list)
    info = result_list(it_result).info;
    arnoldi_info = info.arnoldi_info;
    % The solver accumulates these times over info.iter Arnoldi calls.
    ratios(it_result, :) = [arnoldi_info.t_apply_op, arnoldi_info.t_apply_sketch, ...
        arnoldi_info.t_coeff, arnoldi_info.t_substract, ...
        arnoldi_info.t_update, arnoldi_info.t_other] / info.iter / baseline_time;
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
ylabel(ax, "Average time per Arnoldi call relative to Arnoldi");
title(ax, figure_title, "Interpreter", "none");
legend(ax, time_names, "Interpreter", "none", "Orientation", "horizontal", ...
    "NumColumns", 3, "Location", "southoutside");
exportgraphics(fig, figure_prefix + ".pdf", "ContentType", "vector");
exportgraphics(fig, figure_prefix + ".png", "Resolution", 300);

end
