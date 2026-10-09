function [w, H, h] = arnoldi_last_orth_update(V, m, w, H, h)
%ARNOLDI_LAST_ORTH_UPDATE Orthogonalize the final vector against V.

c = V(:, 1 : m) \ w;
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
