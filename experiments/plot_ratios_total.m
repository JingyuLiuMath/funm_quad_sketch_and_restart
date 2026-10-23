function plot_ratios_total(result_list, figure_prefix, figure_title, baseline_time)
% Normalize timing parts by a shared or per-method baseline total time.

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
group_names = strings(1, numel(result_list));
for it_result = 1 : numel(result_list)
    info = result_list(it_result).info;
    ratios(it_result, :) = [info.t_arnoldi, info.t_funm_init, ...
        info.t_funm_quad, info.t_thick, info.t_update, info.t_other] / baseline_time(it_result);
    bar_labels(it_result) = result_list(it_result).method_name;
    if isfield(result_list, 'group_name')
        group_names(it_result) = result_list(it_result).group_name;
    end
end

y_label = "run time relative to Arn";
plot_ratio_bars(ratios, time_names, bar_labels, figure_prefix, ...
    figure_title, y_label, group_names);

end
