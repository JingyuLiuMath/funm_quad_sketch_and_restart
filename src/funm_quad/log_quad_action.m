function h = log_quad_action(H, rhs, num_quad, nodes, subdiag)
%LOG_QUAD_ACTION One Gauss-Legendre rule for f(z) = log(1+z)/z.
% This is the logarithmic kernel used by the original FUNM_QUAD package.

if num_quad == 1
    t = 0;
    weights = 2;
else
    b_gauss = 0.5 ./ sqrt(1 - (2 * (1:num_quad-1)).^(-2));
    T_gauss = diag(b_gauss, 1) + diag(b_gauss, -1);
    [V_gauss, D_gauss] = eig(T_gauss);
    [t, ind] = sort(diag(D_gauss));
    weights = 2 * V_gauss(1, ind).'.^2;
end
z = -2 ./ (1 + t);
rho = evalnodal(z, nodes, subdiag);

I = eye(size(H));
h = zeros(size(rhs));
for j = 1:num_quad
    h = h + weights(j) * rho(j) * ...
        (((1 + t(j)) * H + 2 * I) \ rhs);
end

end
