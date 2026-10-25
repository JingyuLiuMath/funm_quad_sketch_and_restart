function [v, sv, h, H] = arnoldi_last_shmorth_update(v, h, V, H, sv, SV)

% See the comments in arnoldi_last_hmorth_update.m

m = size(H, 1);
em = zeros(m, 1);
em(m) = 1;
[~, R] = qr([SV, sv], 'econ');
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
