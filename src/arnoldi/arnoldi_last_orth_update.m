function arnoldi_data = arnoldi_last_orth_update(arnoldi_data)

m = size(arnoldi_data.H, 1);
c = arnoldi_data.V(:, 1 : m) \ arnoldi_data.v;
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
