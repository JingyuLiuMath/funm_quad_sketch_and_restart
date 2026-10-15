function arnoldi_data = arnoldi_last_shmorth_update(arnoldi_data)

% See the comments in arnoldi_last_hmorth_update.m

m = size(arnoldi_data.H, 1);
em = zeros(m, 1);
em(m) = 1;
[~, R] = qr([arnoldi_data.SV, arnoldi_data.sv], 'econ');
R11 = R(1 : m, 1 : m);
r = R(1 : m, end);
rho = R(end, end);
M = arnoldi_data.H' * R11' + arnoldi_data.h * (em * r');
c = R11 \ (r + M \ (arnoldi_data.h * abs(rho)^2 * em));

arnoldi_data.v = arnoldi_data.v - arnoldi_data.V(:, 1 : m) * c;
arnoldi_data.sv = arnoldi_data.sv - arnoldi_data.SV * c;
arnoldi_data.H(1 : m, m) = arnoldi_data.H(1 : m, m) + c * arnoldi_data.h;
norm_v = norm(arnoldi_data.v);
if norm_v == 0
    arnoldi_data.v(:) = 0;
    arnoldi_data.sv(:) = 0;
    arnoldi_data.h = 0;
else
    arnoldi_data.v = arnoldi_data.v / norm_v;
    arnoldi_data.sv = arnoldi_data.sv / norm_v;
    arnoldi_data.h = arnoldi_data.h * norm_v;
end

end
