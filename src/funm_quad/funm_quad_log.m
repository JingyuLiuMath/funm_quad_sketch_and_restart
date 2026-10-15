function [f, info] = funm_quad_log(A, b, param)
% f = log(I + A) * A^(-1) \ b

initial_action = @(H, rhs) logm(eye(size(H)) + H) * (H \ rhs);
quad_action = @log_quad_action;
[f, info] = funm_quad_restarted(A, b, param, initial_action, quad_action);

end
