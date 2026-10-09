clear;
close all;
rng(0);

% Two-dimensional Laplacian.
n = 128;
N = n^2;
T = spdiags([-ones(n, 1), 2 * ones(n, 1), -ones(n, 1)], -1:1, n, n);
A = kron(speye(n), T) + kron(T, speye(n));
b = randn(N, 1);
b = b / norm(b);

% Explicit spectrum of A: eigenvectors Q kron Q, eigenvalues d(i) + d(j).
j = 1:n;
Q = sqrt(2 / (n + 1)) * sin((pi / (n + 1)) * (j' * j));
d = 4 * sin(j * pi / (2 * (n + 1))).^2;
lambda = d' + d;
B_hat = Q' * reshape(b, n, n) * Q;
F_exact = Q * (B_hat .* (log1p(lambda) ./ lambda)) * Q';
f_exact = F_exact(:);

param = struct();
param.restart_length = 100;
param.max_restarts = 30;
param.tol = 1e-12;
param.max_num_quad_points = 512;
param.stopping_accuracy = 1e-10;
param.verbose = 0;
param.thick = false;
param.number_thick = 5;
param.exact = f_exact;
param.sketching_mat_type = 'sparse sign';
param.sketching_size_control = 4;
param.cond_tol = 1e3;
param.mode = "adaptive";
param.truncation_length = 2;
param.arnoldi = 'arnoldi';

% Euclidean final-vector correction.
rng(0);
param.update = 'last_orth';
[f_last_orth, info_last_orth] = funm_quad_log(A, b, param);
err_last_orth = norm(f_last_orth - f_exact) / norm(f_exact);
fprintf('last_orth\n');
fprintf('  time: %.1e\n', info_last_orth.t_total);
fprintf('  reason: %s\n', info_last_orth.stop_condition);
fprintf('  iteration number: %d\n', info_last_orth.iter);
fprintf('  relative error: %.1e\n', err_last_orth);
fprintf('  matrix-vector products: %d\n', sum(info_last_orth.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_last_orth.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_last_orth.t_quadrature);

% Harmonic final-vector correction.
rng(0);
param.update = 'last_hmorth';
[f_last_hmorth, info_last_hmorth] = funm_quad_log(A, b, param);
err_last_hmorth = norm(f_last_hmorth - f_exact) / norm(f_exact);
fprintf('last_hmorth\n');
fprintf('  time: %.1e\n', info_last_hmorth.t_total);
fprintf('  reason: %s\n', info_last_hmorth.stop_condition);
fprintf('  iteration number: %d\n', info_last_hmorth.iter);
fprintf('  relative error: %.1e\n', err_last_hmorth);
fprintf('  matrix-vector products: %d\n', sum(info_last_hmorth.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_last_hmorth.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_last_hmorth.t_quadrature);

% Sketched final-vector correction.
rng(0);
param.update = 'last_sorth';
[f_last_sorth, info_last_sorth] = funm_quad_log(A, b, param);
err_last_sorth = norm(f_last_sorth - f_exact) / norm(f_exact);
fprintf('last_sorth\n');
fprintf('  time: %.1e\n', info_last_sorth.t_total);
fprintf('  reason: %s\n', info_last_sorth.stop_condition);
fprintf('  iteration number: %d\n', info_last_sorth.iter);
fprintf('  relative error: %.1e\n', err_last_sorth);
fprintf('  matrix-vector products: %d\n', sum(info_last_sorth.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_last_sorth.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_last_sorth.t_quadrature);

% Sketched harmonic final-vector correction.
rng(0);
param.update = 'last_shmorth';
[f_last_shmorth, info_last_shmorth] = funm_quad_log(A, b, param);
err_last_shmorth = norm(f_last_shmorth - f_exact) / norm(f_exact);
fprintf('last_shmorth\n');
fprintf('  time: %.1e\n', info_last_shmorth.t_total);
fprintf('  reason: %s\n', info_last_shmorth.stop_condition);
fprintf('  iteration number: %d\n', info_last_shmorth.iter);
fprintf('  relative error: %.1e\n', err_last_shmorth);
fprintf('  matrix-vector products: %d\n', sum(info_last_shmorth.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_last_shmorth.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_last_shmorth.t_quadrature);

% Error histories.
figure('Color', 'w');
semilogy(1:info_last_orth.iter, info_last_orth.err_history, '-o', ...
    'LineWidth', 2, 'DisplayName', 'last_orth');
hold on;
semilogy(1:info_last_hmorth.iter, info_last_hmorth.err_history, '-x', ...
    'LineWidth', 2, 'DisplayName', 'last_hmorth');
semilogy(1:info_last_sorth.iter, info_last_sorth.err_history, '-s', ...
    'LineWidth', 2, 'DisplayName', 'last_sorth');
semilogy(1:info_last_shmorth.iter, info_last_shmorth.err_history, '-d', ...
    'LineWidth', 2, 'DisplayName', 'last_shmorth');
legend('Location', 'best', 'Interpreter', 'none');
xlabel('Restart cycle');
ylabel('Relative error');
title('Logarithmic kernel: adaptive');
grid on;
