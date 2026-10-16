function [A, U_dir, lam_dir] = two_dimensional_Laplacian(n)

e = ones(n, 1);
T = spdiags([-e, 2 * e, -e], -1 : 1, n, n);
I = speye(n);
A = kron(I, T) + kron(T, I);

j = (1:n)';
U_dir = sqrt(2 / (n + 1)) * sin((pi / (n + 1)) * (j * j'));
lam_dir = 4 * sin(j * pi / (2 * (n + 1))).^2;

end
