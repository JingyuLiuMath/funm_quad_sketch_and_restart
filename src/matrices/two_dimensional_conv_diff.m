function A = two_dimensional_conv_diff(n, nu)
% Convection-diffusion operator on an n-by-n interior grid.

h = 1 / (n + 1);
e = ones(n, 1);
L = spdiags([-e, 2 * e, -e], -1 : 1, n, n);
C = spdiags([-e, e], [-1, 0], n, n);
G = nu / h^2 * L + C / h;
I = speye(n);
A = kron(I, G) + kron(G.', I);

end
