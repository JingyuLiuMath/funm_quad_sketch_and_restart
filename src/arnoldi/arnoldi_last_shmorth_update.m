function [w, H, h, Sw] = arnoldi_last_shmorth_update(V, m, w, H, h, SV, Sw)
%ARNOLDI_LAST_SHMORTH_UPDATE Use the sketched harmonic projection.

em = zeros(m, 1);
em(m) = 1;
[~, R] = qr([SV, Sw], 'econ');
R11 = R(1 : m, 1 : m);
r = R(1 : m, end);
rho = R(end, end);
M = H' * R11' + h * (em * r');
c = R11 \ (r + M \ (h * rho^2 * em));

w = w - V(:, 1 : m) * c;
Sw = Sw - SV * c;
H(1 : m, m) = H(1 : m, m) + c * h;
norm_w = norm(Sw);
if norm_w == 0
    w(:) = 0;
    h = 0;
else
    w = w / norm_w;
    h = h * norm_w;
end

end
