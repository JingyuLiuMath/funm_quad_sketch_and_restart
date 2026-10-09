function [w, H, h] = arnoldi_last_hmorth_update(V, m, w, H, h)
%ARNOLDI_LAST_HMORTH_UPDATE Make the final vector orthogonal to A*V.
% The Arnoldi relation and a QR factorization avoid forming A*V.

em = zeros(m, 1);
em(m) = 1;
[~, R] = qr([V(:, 1 : m), w], 'econ');
R11 = R(1 : m, 1 : m);
r = R(1 : m, end);
rho = R(end, end);
M = H' * R11' + h * (em * r');
c = R11 \ (r + M \ (h * abs(rho)^2 * em));

w = w - V(:, 1 : m) * c;
H(1 : m, m) = H(1 : m, m) + c * h;
norm_w = norm(w);
if norm_w == 0
    w(:) = 0;
    h = 0;
else
    w = w / norm_w;
    h = h * norm_w;
end

end
