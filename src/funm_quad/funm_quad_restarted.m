function [f, info] = funm_quad_restarted(A, b, ...
    param, initial_action, quad_action)
%FUNM_QUAD_RESTARTED Common restart loop; all work arrays are local.

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
info.update_norm = [];
info.num_quadpoints = [];
info.quad_converged = [];

info.t_arnoldi = 0;
info.t_initial_action = 0;
info.t_quadrature = 0;
info.t_thick = 0;
info.t_update = 0;
info.t_other = 0;
info.t_total = 0;
info.err_history = [];
info.stop_condition = "maximal # of restarts exceeded";

t_total_start = tic;
if norm(b) == 0
    info.stop_condition = "zero right-hand side";
    info.t_total = toc(t_total_start);
    info.t_other = info.t_total;
    return
end

% A fixed sketched method uses this same S for every restart cycle.
S = [];
if param.mode == "fixed" && ~strcmp(param.arnoldi, 'arnoldi')
    S = sketching_mat(param.sketching_size, N, param.sketching_mat_type);
end

m = param.restart_length;
num_quad = min(32, param.max_num_quad_points);

v = b;
beta_acc = 1;
nodes = [];
subdiag = [];

for k = 1 : param.max_restarts
    ell = 0;
    H_init = [];
    if k > 1 && param.thick
        t_thick_start = tic;
        [V_keep, H_init, nodes, ell] = ...
            thick_restarting(V, H, eta, nodes, param.number_thick);
        info.t_thick = info.t_thick + toc(t_thick_start);
    end
    t_arnoldi_start = tic;
    if k == 1
        V = zeros(N, m);
    else
        V = zeros(N, m + ell);
    end
    if ell > 0
        V(:, 1 : ell) = V_keep;
    end
    V(:, ell + 1) = v;
    if param.mode == "adaptive"
        [V, v, H, eta, beta, succ_flag, num_oracle, S] = ...
            arnoldi_adaptive(S, A, V, m + ell, H_init, ell + 1, param);
        m_now = size(H, 1);
        if k == 1
            m = m_now;
        end
    else
        switch param.arnoldi
            case 'arnoldi'
                [V, v, H, eta, beta, succ_flag, num_oracle] = ...
                    arnoldi_fix(A, V, m + ell, H_init, ell + 1, param);
            case 'sarnoldi'
                [V, v, H, eta, beta, succ_flag, num_oracle] = ...
                    sarnoldi_fix(S, A, V, m + ell, H_init, ell + 1, param);
            case 'shmarnoldi'
                [V, v, H, eta, beta, succ_flag, num_oracle] = ...
                    shmarnoldi_fix(S, A, V, m + ell, H_init, ell + 1, param);
        end
        m_now = size(H, 1);
    end
    info.t_arnoldi = info.t_arnoldi + toc(t_arnoldi_start);
    info.iter = k;
    info.dim(k) = m_now;
    info.num_thick(k) = ell;
    info.num_oracle(k) = num_oracle;
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

    beta_acc = beta_acc * beta;
    rhs = zeros(m_now, 1);
    rhs(ell + 1) = beta_acc;
    if k == 1
        t_initial_action_start = tic;
        h = initial_action(H, rhs);
        info.t_initial_action = info.t_initial_action + toc(t_initial_action_start);
        num_quad_used = 0;
        quad_converged = true;
    else
        t_quadrature_start = tic;
        [h, num_quad_used, quad_converged, num_quad] = ...
            refine_quadrature(H, rhs, nodes, subdiag, num_quad, param, quad_action);
        info.t_quadrature = info.t_quadrature + toc(t_quadrature_start);
    end

    t_update_start = tic;
    f_update = V * h;
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
        info.stop_condition = "relative update norm below stopping accuracy";
        break
    end
    if k < param.max_restarts
        nodes = [nodes(1 : end - ell); eig(H)];
        H_new = H(ell + 1 : m_now, ell + 1 : m_now);
        subdiag = [subdiag; H_new(2 : size(H_new, 1) + 1 : end).'; eta];
    end
end

info.t_total = toc(t_total_start);
info.t_other = info.t_total - (info.t_arnoldi + ...
    info.t_initial_action + info.t_quadrature + info.t_thick + info.t_update);

end
