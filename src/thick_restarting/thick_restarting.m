function [V, H, SV, nodes, num_thick] = ...
    thick_restarting(h, V, H, SV, nodes, number_thick, use_sketch)
%THICK_RESTARTING Retain a Schur subspace and seed the next Arnoldi cycle.
% For Arnoldi and sketched Arnoldi; harmonic duals are not supported.
% H is returned as a (num_thick + 1)-by-num_thick seed matrix.
% The caller keeps the remainder v and its optional sketch sv unchanged.

m_prev = size(H, 1);
if isreal(H)
    [U, T] = schur(H, 'real');
else
    [U, T] = schur(H, 'complex');
end
[num_thick, U, ~, evT] = thick_quad(U, T, ordeig(T), number_thick);
U = U(:, 1 : num_thick);
V = V * U;
H = [U' * H * U; h * U(end, :)];
if use_sketch
    SV = SV * U;
end

% Retained roots move to the tail and are replaced by the next cycle's
% projected roots when its quadrature history is appended.
nodes = [nodes(1 : end - m_prev); evT(num_thick + 1 : end); evT(1 : num_thick)];

end
