function [v, sv, h, H] = arnoldi_last_hmorth_update(v, h, V, H, sv, SV)

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

m = size(H, 1);
em = zeros(m, 1);
em(m) = 1;
[~, R] = qr([V, v], 'econ');
R11 = R(1 : m, 1 : m);
r = R(1 : m, end);
rho = R(end, end);
M = H' * R11' + h * (em * r');
c = R11 \ (r + M \ (h * abs(rho)^2 * em));

v = v - V * c;
sv = sv - SV * c;
H(1 : m, m) = H(1 : m, m) + c * h;
norm_v = norm(v);
if norm_v == 0
    v(:) = 0;
    sv(:) = 0;
    h = 0;
else
    v = v / norm_v;
    sv = sv / norm_v;
    h = h * norm_v;
end

end
