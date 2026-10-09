function [f, info] = funm_quad_invsqrt(A, b, param)
%FUNM_QUAD_INVSQRT Approximate A^(-1/2)*b by restarted Krylov quadrature.

initial_action = @(H, rhs) sqrtm(H) \ rhs;
quad_action = @invsqrt_quad_action;
[f, info] = funm_quad_restarted(A, b, param, initial_action, quad_action);

end
