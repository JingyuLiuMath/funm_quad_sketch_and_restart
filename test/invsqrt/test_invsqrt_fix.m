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
param.sketch_dim = 4 * param.krylov_dim;
param.mode = "fixed";
param.truncation_length = inf;
param.update = "";

% Standard Arnoldi.
rng(0);
param.arnoldi = 'arnoldi';
[f_arnoldi, info_arnoldi] = funm_quad_invsqrt(A, b, param);
fprintf('arnoldi\n');
fprintf('  time: %.1e\n', info_arnoldi.t_total);
fprintf('  reason: %s\n', info_arnoldi.stop_condition);
fprintf('  iteration number: %d\n', info_arnoldi.iter);
fprintf('  relative error: %.1e\n', info_arnoldi.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_arnoldi.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_arnoldi.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_arnoldi.t_funm_quad);

% Harmonic Arnoldi.
rng(0);
param.arnoldi = 'hmarnoldi';
[f_hmarnoldi, info_hmarnoldi] = funm_quad_invsqrt(A, b, param);
fprintf('hmarnoldi\n');
fprintf('  time: %.1e\n', info_hmarnoldi.t_total);
fprintf('  reason: %s\n', info_hmarnoldi.stop_condition);
fprintf('  iteration number: %d\n', info_hmarnoldi.iter);
fprintf('  relative error: %.1e\n', info_hmarnoldi.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_hmarnoldi.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_hmarnoldi.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_hmarnoldi.t_funm_quad);

% Sketched Arnoldi.
rng(0);
param.arnoldi = 'sarnoldi';
[f_sarnoldi, info_sarnoldi] = funm_quad_invsqrt(A, b, param);
fprintf('sarnoldi\n');
fprintf('  time: %.1e\n', info_sarnoldi.t_total);
fprintf('  reason: %s\n', info_sarnoldi.stop_condition);
fprintf('  iteration number: %d\n', info_sarnoldi.iter);
fprintf('  relative error: %.1e\n', info_sarnoldi.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_sarnoldi.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_sarnoldi.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_sarnoldi.t_funm_quad);

% Sketched harmonic Arnoldi.
rng(0);
param.arnoldi = 'shmarnoldi';
[f_shmarnoldi, info_shmarnoldi] = funm_quad_invsqrt(A, b, param);
fprintf('shmarnoldi\n');
fprintf('  time: %.1e\n', info_shmarnoldi.t_total);
fprintf('  reason: %s\n', info_shmarnoldi.stop_condition);
fprintf('  iteration number: %d\n', info_shmarnoldi.iter);
fprintf('  relative error: %.1e\n', info_shmarnoldi.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_shmarnoldi.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_shmarnoldi.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_shmarnoldi.t_funm_quad);

% Error histories.
figure('Color', 'w');
semilogy(1:info_arnoldi.iter, info_arnoldi.err_history, '-o', ...
    'LineWidth', 2, 'DisplayName', 'arnoldi');
hold on;
semilogy(1:info_hmarnoldi.iter, info_hmarnoldi.err_history, '-x', ...
    'LineWidth', 2, 'DisplayName', 'hmarnoldi');
semilogy(1:info_sarnoldi.iter, info_sarnoldi.err_history, '-s', ...
    'LineWidth', 2, 'DisplayName', 'sarnoldi');
semilogy(1:info_shmarnoldi.iter, info_shmarnoldi.err_history, '-d', ...
    'LineWidth', 2, 'DisplayName', 'shmarnoldi');
legend('Location', 'best', 'Interpreter', 'none');
xlabel('Restart cycle');
ylabel('Relative error');
title('Inverse square root: fixed');
