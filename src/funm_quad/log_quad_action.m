function h = log_quad_action(H, rhs, num_quad, nodes, subdiag)

% Quadrature for the restarted log(1+x)/x action.
%
% The principal branch has the Stieltjes representation
%   f(x) = log(1+x)/x = integral_{1}^{inf} 1/(s*(x+s)) ds,
% with f(0) = 1 by continuity.  Define the accumulated restart factor
%   rho(z) = prod_{ell} (subdiag(ell)/(z-nodes(ell))).
% This routine approximates the correction integral
%   h = integral_{1}^{inf} rho(-s)*(H+s*I)^(-1)*rhs/s ds.
% The spectrum of H and the nodal poles must lie outside (-inf, -1].
% If nodes is empty, rho = 1 and the integral equals f(H)*rhs.
%
% Set s = 2/(1+t), z(t) = -s, and B(t) = (1+t)*H + 2*I.
% The transformed integral is
%   h = integral_{-1}^{1} rho(z(t))*B(t)^(-1)*rhs dt.
% With n = num_quad, Gauss-Legendre quadrature uses nodes t_j and weights
% w_j for the unit weight on [-1,1], giving
%   h approx sum_{j=1}^{n} w_j*rho(z(t_j))*B(t_j)^(-1)*rhs.

if num_quad == 1
    t = 0;
    weights = 2;
else
    % Golub-Welsch construction of Gauss-Legendre nodes and weights.
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
