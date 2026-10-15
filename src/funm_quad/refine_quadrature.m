function [h, num_quad_used, converged, num_quad_next] = refine_quadrature( ...
    H, rhs, nodes, subdiag, num_quad, param, quad_action)

num_quad = min(max(num_quad, 16), param.max_num_quad_points);
% The half-size coarse rule is used only to estimate the quadrature error.
num_quad_coarse = max(8, floor(num_quad / 2));
h_coarse = quad_action(H, rhs, num_quad_coarse, nodes, subdiag);
num_quad_used = num_quad;
while true
    h = quad_action(H, rhs, num_quad_used, nodes, subdiag);
    delta = norm(h - h_coarse);
    hnorm = norm(h);
    converged = delta <= param.quad_tol * hnorm ;
    if converged || num_quad_used >= param.max_num_quad_points
        break
    end
    h_coarse = h;
    num_quad_coarse = num_quad_used;
    num_quad_used = min(2 * num_quad_used, param.max_num_quad_points);
end

num_quad_next = max(16, num_quad_coarse);

end
