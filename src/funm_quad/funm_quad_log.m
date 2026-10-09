function [f, info] = funm_quad_log(A, b, param)
%FUNM_QUAD_LOG Approximate g(A)*b, where g(z)=log(1+z)/z.

initial_action = @(H, rhs) logm(eye(size(H)) + H) * (H \ rhs);
quad_action = @log_quad_action;
[f, info] = funm_quad_restarted(A, b, param, initial_action, quad_action);

end
