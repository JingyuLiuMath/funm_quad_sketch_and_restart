clear;
close all;
rng(0);

% Two-dimensional Laplacian.
n = 128;
N = n^2;
[A, U_dir, lam_dir] = two_dimensional_Laplacian(n);
b = randn(N, 1);
b = b / norm(b);

% Reference solution for g(z) = log(1 + z) / z.
lam = lam_dir + lam_dir';
B_proj = U_dir' * reshape(b, n, n) * U_dir;
F_exact = U_dir * (B_proj .* (log1p(lam) ./ lam)) * U_dir';
f_exact = F_exact(:);

param = struct();
param.krylov_dim = 100;
param.max_restarts = 30;
param.quad_tol = 1e-12;
param.max_num_quad_points = 512;
param.stopping_accuracy = 1e-10;
param.number_thick = 5;
param.exact = f_exact;
param.sketching_size_control = 4;
param.cond_tol = 1e3;
param.mode = "adaptive";
param.truncation_length = 1;
param.arnoldi = 'arnoldi';

% last_orth.
rng(0);
param.thick = false;
param.update = 'last_orth';
[f_last_orth, info_last_orth] = funm_quad_log(A, b, param);
fprintf('last_orth\n');
fprintf('  time: %.1e\n', info_last_orth.t_total);
fprintf('  reason: %s\n', info_last_orth.stop_condition);
fprintf('  iteration number: %d\n', info_last_orth.iter);
fprintf('  relative error: %.1e\n', info_last_orth.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_last_orth.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_last_orth.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_last_orth.t_funm_quad);
fprintf('  thick restart time: %.1e\n', info_last_orth.t_thick);

% last_sorth.
rng(0);
param.thick = false;
param.update = 'last_sorth';
[f_last_sorth, info_last_sorth] = funm_quad_log(A, b, param);
fprintf('last_sorth\n');
fprintf('  time: %.1e\n', info_last_sorth.t_total);
fprintf('  reason: %s\n', info_last_sorth.stop_condition);
fprintf('  iteration number: %d\n', info_last_sorth.iter);
fprintf('  relative error: %.1e\n', info_last_sorth.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_last_sorth.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_last_sorth.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_last_sorth.t_funm_quad);
fprintf('  thick restart time: %.1e\n', info_last_sorth.t_thick);

% last_orth (thick restart).
rng(0);
param.thick = true;
param.update = 'last_orth';
[f_last_orth_thick, info_last_orth_thick] = funm_quad_log(A, b, param);
fprintf('last_orth (thick restart)\n');
fprintf('  time: %.1e\n', info_last_orth_thick.t_total);
fprintf('  reason: %s\n', info_last_orth_thick.stop_condition);
fprintf('  iteration number: %d\n', info_last_orth_thick.iter);
fprintf('  relative error: %.1e\n', info_last_orth_thick.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_last_orth_thick.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_last_orth_thick.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_last_orth_thick.t_funm_quad);
fprintf('  thick restart time: %.1e\n', info_last_orth_thick.t_thick);

% last_sorth (thick restart).
rng(0);
param.thick = true;
param.update = 'last_sorth';
[f_last_sorth_thick, info_last_sorth_thick] = funm_quad_log(A, b, param);
fprintf('last_sorth (thick restart)\n');
fprintf('  time: %.1e\n', info_last_sorth_thick.t_total);
fprintf('  reason: %s\n', info_last_sorth_thick.stop_condition);
fprintf('  iteration number: %d\n', info_last_sorth_thick.iter);
fprintf('  relative error: %.1e\n', info_last_sorth_thick.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_last_sorth_thick.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_last_sorth_thick.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_last_sorth_thick.t_funm_quad);
fprintf('  thick restart time: %.1e\n', info_last_sorth_thick.t_thick);

% Error histories.
figure('Color', 'w');
semilogy(1:info_last_orth.iter, info_last_orth.err_history, '-o', ...
    'Color', [0, 0.4470, 0.7410], 'LineWidth', 2, ...
    'DisplayName', 'last_orth');
hold on;
semilogy(1:info_last_sorth.iter, info_last_sorth.err_history, '-s', ...
    'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, ...
    'DisplayName', 'last_sorth');
semilogy(1:info_last_orth_thick.iter, info_last_orth_thick.err_history, '--o', ...
    'Color', [0, 0.4470, 0.7410], 'LineWidth', 2, ...
    'DisplayName', 'last_orth (thick restart)');
semilogy(1:info_last_sorth_thick.iter, info_last_sorth_thick.err_history, '--s', ...
    'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, ...
    'DisplayName', 'last_sorth (thick restart)');
legend('Location', 'best', 'Interpreter', 'none');
xlabel('Restart cycle');
ylabel('Relative error');
title('Logarithmic kernel: adaptive with/without thick restart');

