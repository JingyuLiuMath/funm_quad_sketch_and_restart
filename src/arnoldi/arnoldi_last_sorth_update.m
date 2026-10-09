function [w, H, h, Sw] = arnoldi_last_sorth_update(V, m, w, H, h, SV, Sw)
%ARNOLDI_LAST_SORTH_UPDATE Orthogonalize using the supplied sketch images.

c = SV \ Sw;
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
