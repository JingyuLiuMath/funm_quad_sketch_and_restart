function [f, info] = funm_quad_restarted(A, b, ...
    param, initial_action, quad_action)

N = size(b, 1);

has_exact = isfield(param, 'exact') && ~isempty(param.exact);
if has_exact
    norm_exact = norm(param.exact);
end

f = zeros(size(b));

info = struct();
info.param = param;
info.iter = 0;
info.dim = [];
info.num_thick = [];
info.num_oracle = [];
info.num_resketch = [];
info.update_norm = [];
info.num_quadpoints = [];
info.quad_converged = [];

info.t_arnoldi = 0;
info.t_funm_init = 0;
info.t_funm_quad = 0;
info.t_thick = 0;
info.t_update = 0;
info.t_other = 0;
info.t_total = 0;
info.err_history = [];
info.stop_condition = "maximal number of restarts exceeded";

t_total_start = tic;
if norm(b) == 0
    info.stop_condition = "zero right-hand side";
    info.t_total = toc(t_total_start);
    info.t_other = info.t_total;
    return
end

is_adaptive = param.mode == "adaptive";
use_sketch = is_adaptive || param.arnoldi == "sarnoldi" ...
    || param.arnoldi == "shmarnoldi";
use_harmonic = ~is_adaptive && ...
    (param.arnoldi == "hmarnoldi" || param.arnoldi == "shmarnoldi");

S_op = [];
if is_adaptive
    S_op = sparse_sign(30, N);
elseif use_sketch
    S_op = sparse_sign(param.sketch_dim, N);
end

m = param.krylov_dim;
num_quad = min(32, param.max_num_quad_points);

arnoldi_data = struct();
arnoldi_data.v = b;
if use_sketch
    arnoldi_data.sv = S_op(b);
end
beta_acc = 1;
nodes = [];
subdiag = [];

for k = 1 : param.max_restarts
    num_thick = 0;
    if k > 1 && param.thick
        t_thick_start = tic;
        [arnoldi_data, nodes, num_thick] = ...
            thick_restarting(arnoldi_data, nodes, param.number_thick, use_sketch);
        info.t_thick = info.t_thick + toc(t_thick_start);
    end

    t_arnoldi_start = tic;
    m_total = m + num_thick;
    start_ind = num_thick + 1;
    if num_thick == 0
        arnoldi_data.V = zeros(N, 0);
        arnoldi_data.H = [];
        if use_harmonic
            arnoldi_data.Z = zeros(N, 0);
        end
        if use_sketch
            s = size(arnoldi_data.sv, 1);
            arnoldi_data.SV = zeros(s, 0);
            if use_harmonic
                arnoldi_data.SZ = zeros(s, 0);
            end
        end
    end

    num_resketch = 0;
    if is_adaptive
        [S_op, m_now, arnoldi_data, succ_flag, num_oracle, num_resketch] = ...
            arnoldi_adaptive(S_op, A, arnoldi_data, m_total, start_ind, param);
        if k == 1
            m = m_now;
        end
    else
        switch param.arnoldi
            case "arnoldi"
                [arnoldi_data, succ_flag, num_oracle] = ...
                    arnoldi_fix(A, arnoldi_data, m_total, start_ind, param);
            case "hmarnoldi"
                [arnoldi_data, succ_flag, num_oracle] = ...
                    hmarnoldi_fix(A, arnoldi_data, m_total, start_ind, param);
            case "sarnoldi"
                [arnoldi_data, succ_flag, num_oracle] = ...
                    sarnoldi_fix(S_op, A, arnoldi_data, m_total, start_ind, param);
            case "shmarnoldi"
                [arnoldi_data, succ_flag, num_oracle] = ...
                    shmarnoldi_fix(S_op, A, arnoldi_data, m_total, start_ind, param);
        end
        m_now = size(arnoldi_data.H, 1);
    end
    info.t_arnoldi = info.t_arnoldi + toc(t_arnoldi_start);
    info.iter = k;
    info.dim(k) = m_now;
    info.num_thick(k) = num_thick;
    info.num_oracle(k) = num_oracle;
    info.num_resketch(k) = num_resketch;
    if ~succ_flag
        info.update_norm(k) = 0;
        info.num_quadpoints(k) = 0;
        info.quad_converged(k) = false;
        info.stop_condition = "Arnoldi failed";
        if has_exact
            info.err_history(k) = norm(f - param.exact) / norm_exact;
        end
        break
    end

    beta_acc = beta_acc * arnoldi_data.beta;
    rhs = zeros(m_now, 1);
    rhs(start_ind) = beta_acc;
    if k == 1
        t_funm_init_start = tic;
        f_small = initial_action(arnoldi_data.H, rhs);
        info.t_funm_init = info.t_funm_init + toc(t_funm_init_start);
        num_quad_used = 0;
        quad_converged = true;
    else
        t_funm_quad_start = tic;
        [f_small, num_quad_used, quad_converged, num_quad] = ...
            refine_quadrature(arnoldi_data.H, rhs, nodes, subdiag, ...
            num_quad, param, quad_action);
        info.t_funm_quad = info.t_funm_quad + toc(t_funm_quad_start);
    end

    t_update_start = tic;
    f_update = arnoldi_data.V * f_small;
    f = f + f_update;
    info.update_norm(k) = norm(f_update);
    conv_flag = info.update_norm(k) / norm(f) < param.stopping_accuracy;
    info.t_update = info.t_update + toc(t_update_start);
    info.num_quadpoints(k) = num_quad_used;
    info.quad_converged(k) = quad_converged;
    if has_exact
        info.err_history(k) = norm(f - param.exact) / norm_exact;
    end
    if conv_flag
        info.stop_condition = "converged";
        break
    end
    if k < param.max_restarts
        % Retained roots are replaced; only new-step factors are appended.
        nodes = [nodes(1 : end - num_thick); eig(arnoldi_data.H)];
        H_new = arnoldi_data.H(start_ind : m_now, start_ind : m_now);
        num_new = m_now - num_thick;
        subdiag = [subdiag; ...
            H_new(2 : (num_new + 1) : end).'; arnoldi_data.h];
    end
end

info.t_total = toc(t_total_start);
info.t_other = info.t_total - (info.t_arnoldi + ...
    info.t_funm_init + info.t_funm_quad + info.t_thick + info.t_update);

end
