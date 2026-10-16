clear;
close all;
rng(0);

% Two-dimensional Laplacian.
n = 128;
N = n^2;
[A, U_dir, lam_dir] = two_dimensional_Laplacian(n);
b = randn(N, 1);
b = b / norm(b);

% Reference solution from the one-dimensional sine basis and spectrum.
lam = lam_dir + lam_dir';
B_proj = U_dir' * reshape(b, n, n) * U_dir;
F_exact = U_dir * (B_proj ./ sqrt(lam)) * U_dir';
f_exact = F_exact(:);

param = struct();
param.krylov_dim = 100;
param.max_restarts = 30;
param.quad_tol = 1e-12;
param.max_num_quad_points = 512;
param.stopping_accuracy = 1e-10;
param.thick = false;
param.number_thick = 5;
param.exact = f_exact;
param.sketching_size_control = 4;
param.mode = "adaptive";
param.arnoldi = 'arnoldi';
param.truncation_length = 1;
param.cond_tol = 1e3;

% Standard last-vector correction.
rng(0);
param.update = 'last_orth';
[f_last_orth, info_last_orth] = funm_quad_invsqrt(A, b, param);
fprintf('last_orth\n');
fprintf('  time: %.1e\n', info_last_orth.t_total);
fprintf('  reason: %s\n', info_last_orth.stop_condition);
fprintf('  iteration number: %d\n', info_last_orth.iter);
fprintf('  relative error: %.1e\n', info_last_orth.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_last_orth.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_last_orth.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_last_orth.t_funm_quad);

% Harmonic last-vector correction.
rng(0);
param.update = 'last_hmorth';
[f_last_hmorth, info_last_hmorth] = funm_quad_invsqrt(A, b, param);
fprintf('last_hmorth\n');
fprintf('  time: %.1e\n', info_last_hmorth.t_total);
fprintf('  reason: %s\n', info_last_hmorth.stop_condition);
fprintf('  iteration number: %d\n', info_last_hmorth.iter);
fprintf('  relative error: %.1e\n', info_last_hmorth.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_last_hmorth.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_last_hmorth.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_last_hmorth.t_funm_quad);

% Sketched last-vector correction.
rng(0);
param.update = 'last_sorth';
[f_last_sorth, info_last_sorth] = funm_quad_invsqrt(A, b, param);
fprintf('last_sorth\n');
fprintf('  time: %.1e\n', info_last_sorth.t_total);
fprintf('  reason: %s\n', info_last_sorth.stop_condition);
fprintf('  iteration number: %d\n', info_last_sorth.iter);
fprintf('  relative error: %.1e\n', info_last_sorth.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_last_sorth.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_last_sorth.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_last_sorth.t_funm_quad);

% Sketched harmonic last-vector correction.
rng(0);
param.update = 'last_shmorth';
[f_last_shmorth, info_last_shmorth] = funm_quad_invsqrt(A, b, param);
fprintf('last_shmorth\n');
fprintf('  time: %.1e\n', info_last_shmorth.t_total);
fprintf('  reason: %s\n', info_last_shmorth.stop_condition);
fprintf('  iteration number: %d\n', info_last_shmorth.iter);
fprintf('  relative error: %.1e\n', info_last_shmorth.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_last_shmorth.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_last_shmorth.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_last_shmorth.t_funm_quad);

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
title('Inverse square root: adaptive');
