function [arnoldi_data, nodes, num_thick] = ...
    thick_restarting(arnoldi_data, nodes, number_thick, use_sketch)
%THICK_RESTARTING Retain a Schur subspace and seed the next Arnoldi cycle.
% For Arnoldi and sketched Arnoldi; harmonic duals are not supported.
% H is returned as a (num_thick + 1)-by-num_thick seed matrix.
% The remainder v and its optional sketch sv retain their incoming scaling.

m_prev = size(arnoldi_data.H, 1);
if isreal(arnoldi_data.H)
    [U, T] = schur(arnoldi_data.H, 'real');
else
    [U, T] = schur(arnoldi_data.H, 'complex');
end
[num_thick, U, ~, evT] = thick_quad(U, T, ordeig(T), number_thick);
U = U(:, 1 : num_thick);
arnoldi_data.V = arnoldi_data.V * U;
arnoldi_data.H = [U' * arnoldi_data.H * U; arnoldi_data.h * U(end, :)];
if use_sketch
    arnoldi_data.SV = arnoldi_data.SV * U;
end

% Retained roots move to the tail and are replaced by the next cycle's
% projected roots when its quadrature history is appended.
nodes = [nodes(1 : end - m_prev); evT(num_thick + 1 : end); evT(1 : num_thick)];

end
