function f = two_dimensional_conv_diff_log(n, nu, b)
% Compute log(I + A) * b for A = two_dimensional_conv_diff(n, nu).
% Assumes nu > 0 and an n^2-entry vector b; uses only n-by-n arrays.
%
% The convection-diffusion matrix has the Kronecker sum structure
%   A = kron(I, G) + kron(G.', I).
%
% Under MATLAB column-major vectorization, x = vec(X) satisfies
%   A * x = vec(G * X + X * G).
% Thus, the action of A can be represented on an n-by-n array.
%
% Suppose the one-dimensional matrix has the eigendecomposition
%   G = V * Lambda * V^(-1),    Lambda = diag(lambda_1, ..., lambda_n).
% With X_hat = V^(-1) * X * V,
%   G * X + X * G = V * (Lambda * X_hat + X_hat * Lambda) * V^(-1).
% Each entry (i, j) of X_hat is multiplied by lambda_i + lambda_j.
% The matrix logarithm therefore acts by elementwise multiplication:
%   f = vec(V * (H .* (V^(-1) * B * V)) * V^(-1)),
%   H_ij = log(1 + lambda_i + lambda_j).
% Here B = reshape(b, n, n), and .* denotes elementwise multiplication.
% For b = ones(n^2, 1) / n, this gives B = ones(n) / n.
%
% The eigendecomposition of G is explicit. Let
%   a = nu / h^2,    c = 1 / h,    r = sqrt((a + c) / a),
%   D = diag(1, r, ..., r^(n - 1)).
% Then S = D^(-1) * G * D is a symmetric tridiagonal Toeplitz matrix:
%   S = tridiag(-sqrt(a * (a + c)), 2 * a + c, -sqrt(a * (a + c))).
% Its orthonormal eigenvectors form the sine matrix
%   Q_jk = sqrt(2 / (n + 1)) * sin(pi * j * k / (n + 1)).
% Hence G = D * Q * Lambda * Q.' * D^(-1), with V = D * Q.
%
% The code centers the diagonal scaling D without changing this similarity
% transform.

h = 1 / (n + 1);
a = nu / h^2;
c = 1 / h;
k = (1 : n).';

% Center the diagonal scaling without changing the similarity transform.
gamma = 0.5 * log1p(c / a);
d = exp((k - (n + 1) / 2) * gamma);
% Multiplication by Q can be implemented using a fast orthonormal DST-I.
Q = sqrt(2 / (n + 1)) * sin((pi / (n + 1)) * (k * k.'));

% G = D * Q * diag(lambda) * Q.' / D; avoid cancellation in lambda.
t = sqrt(a) * sqrt(a + c);
lambda = (c / (sqrt(a + c) + sqrt(a)))^2 ...
    + 4 * t * sin(pi * k / (2 * (n + 1))).^2;

% A * vec(B) = vec(G * B + B * G), so transform B on both sides.
B = reshape(b, n, n);
C = Q.' * ((B ./ d) .* d.') * Q;
W = log1p(lambda + lambda.') .* C;

% Transform back to the original coordinates and MATLAB column ordering.
F = (d .* (Q * W * Q.')) ./ d.';
f = F(:);

end
