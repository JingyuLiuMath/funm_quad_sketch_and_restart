function plot_ratios_arnoldi(result_list, figure_prefix, figure_title, baseline_time)
% Compare fixed and adaptive methods using average time per Arnoldi call.
% Baselines are shared or supplied separately for each method.

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

y_label = "run time relative to Arn";
plot_ratio_bars(ratios, time_names, bar_labels, figure_prefix, ...
    figure_title, y_label, group_names);

end
