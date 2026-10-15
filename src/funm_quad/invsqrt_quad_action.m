function h = invsqrt_quad_action(H, rhs, num_quad, nodes, subdiag)

% Quadrature for the restarted inverse-square-root action.
%
% The principal branch has the Stieltjes representation
%   f(x) = x^(-1/2) = (1/pi) * integral_{0}^{inf} 1/(sqrt(s)*(x+s)) ds.
% Define the accumulated restart factor
%   rho(z) = prod_{ell} (subdiag(ell)/(z-nodes(ell))).
% This routine approximates the correction integral
%   h = (1/pi) * integral_{0}^{inf} rho(-s)*(H+s*I)^(-1)*rhs/sqrt(s) ds.
% The spectrum of H and the nodal poles must lie outside (-inf, 0].
% If nodes is empty, rho = 1 and the integral equals H^(-1/2)*rhs.
%
% Set s = (1-t)/(1+t), z(t) = -s, and B(t) = (1-t)*I + (1+t)*H.
% The transformed integral is
%   h = (2/pi) * integral_{-1}^{1} rho(z(t))*B(t)^(-1)*rhs/sqrt(1-t^2) dt.
% With n = num_quad, Gauss-Chebyshev quadrature of the first kind uses
%   t_j = cos((2*j-1)*pi/(2*n)),  w_j = pi/n,  j = 1,...,n,
% for the weight 1/sqrt(1-t^2), giving
%   h approx (2/pi) * sum_{j=1}^{n} w_j*rho(z(t_j))*B(t_j)^(-1)*rhs.

t = cos((2 * (1:num_quad).' - 1) * pi / (2 * num_quad));
weights = (pi / num_quad) * ones(num_quad, 1);
z = -(1 - t) ./ (1 + t);
rho = evalnodal(z, nodes, subdiag);

I = eye(size(H));
h = zeros(size(rhs));
% Each solve uses -B(t_j), so the final factor -2/pi gives the sum above.
for j = 1:num_quad
    h = h + weights(j) * rho(j) * ...
        ((-(1 - t(j)) * I - (1 + t(j)) * H) \ rhs);
end
h = -2 / pi * h;

end
