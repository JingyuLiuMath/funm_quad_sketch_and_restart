function h = invsqrt_quad_action(H, rhs, num_quad, nodes, subdiag)
%INVSQRT_QUAD_ACTION One Gauss-Jacobi rule for the inverse-square-root restart.
% The variable transformation uses unit scale.

t = cos((2 * (1:num_quad).' - 1) * pi / (2 * num_quad));
weights = (pi / num_quad) * ones(num_quad, 1);
z = -(1 - t) ./ (1 + t);
rho = evalnodal(z, nodes, subdiag);

I = eye(size(H));
h = zeros(size(rhs));
for j = 1:num_quad
    h = h + weights(j) * rho(j) * ...
        ((-(1 - t(j)) * I - (1 + t(j)) * H) \ rhs);
end
h = -2 / pi * h;

end
