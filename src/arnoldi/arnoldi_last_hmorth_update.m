function arnoldi_data = arnoldi_last_hmorth_update(arnoldi_data)

% Let V, v, H, h denote the input data, with real h and
%     A*V = V*H + h*v*em'.
% The harmonic correction chooses c so that (A*V)'*(v - V*c) = 0:
%     ((A*V)'*V)*c = (A*V)'*v.
%
% Factor [V, v] = [Q1, q]*[R11, r; 0, rho] and define
%     B = R11*H + h*r*em'.
% The Arnoldi relation then gives
%     A*V = Q1*B + q*(h*rho*em'),
%     M = B' = H'*R11' + h*em*r',
%     (A*V)'*V = M*R11,
%     (A*V)'*v = M*r + h*abs(rho)^2*em.
% For nonsingular R11 and M, the correction is therefore
%     c = R11\(r + M\(h*abs(rho)^2*em)).
%

m = size(arnoldi_data.H, 1);
em = zeros(m, 1);
em(m) = 1;
[~, R] = qr([arnoldi_data.V(:, 1 : m), arnoldi_data.v], 'econ');
R11 = R(1 : m, 1 : m);
r = R(1 : m, end);
rho = R(end, end);
M = arnoldi_data.H' * R11' + arnoldi_data.h * (em * r');
c = R11 \ (r + M \ (arnoldi_data.h * abs(rho)^2 * em));

arnoldi_data.v = arnoldi_data.v - arnoldi_data.V(:, 1 : m) * c;
arnoldi_data.H(1 : m, m) = arnoldi_data.H(1 : m, m) + c * arnoldi_data.h;
norm_v = norm(arnoldi_data.v);
if norm_v == 0
    arnoldi_data.v(:) = 0;
    arnoldi_data.h = 0;
else
    arnoldi_data.v = arnoldi_data.v / norm_v;
    arnoldi_data.h = arnoldi_data.h * norm_v;
end

end
