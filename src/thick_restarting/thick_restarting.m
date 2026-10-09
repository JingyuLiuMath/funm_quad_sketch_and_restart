function [V_keep, H_init, nodes, num_thick] = ...
    thick_restarting(V, H, eta, nodes, number_thick)
%THICK_RESTARTING Prepare the retained basis and seeded Arnoldi relation.
% The current remainder vector is kept by the caller. The last row of
% H_init contains its coefficients for the retained Schur subspace.

m_prev = size(H, 1);
if isreal(H)
    [U, T] = schur(H, 'real');
else
    [U, T] = schur(H, 'complex');
end
[num_thick, U, ~, D] = thick_quad(U, T, ordeig(T), number_thick);
U = U(:, 1 : num_thick);
V_keep = V * U;
H_init = [U' * H * U; eta * U(end, :)];

% Retained roots move to the tail and are replaced by the next cycle's
% projected roots when its quadrature history is appended.
nodes = [nodes(1 : end - m_prev); D(num_thick + 1 : end); D(1 : num_thick)];

end
