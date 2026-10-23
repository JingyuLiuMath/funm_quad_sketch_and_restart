function info = read_median_run(selection, mode, method, truncation_length)
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
