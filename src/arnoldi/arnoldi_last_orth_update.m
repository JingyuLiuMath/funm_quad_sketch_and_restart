function [v, sv, h, H] = arnoldi_last_orth_update(v, h, V, H, sv, SV)

m = size(H, 1);
c = V \ v;
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
