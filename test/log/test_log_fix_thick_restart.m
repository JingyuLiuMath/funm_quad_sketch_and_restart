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
param.sketch_dim = 4 * param.krylov_dim;
param.mode = "fixed";
param.truncation_length = inf;
param.update = "";

% arnoldi.
rng(0);
param.thick = false;
param.arnoldi = 'arnoldi';
[f_arnoldi, info_arnoldi] = funm_quad_log(A, b, param);
fprintf('arnoldi\n');
fprintf('  time: %.1e\n', info_arnoldi.t_total);
fprintf('  reason: %s\n', info_arnoldi.stop_condition);
fprintf('  iteration number: %d\n', info_arnoldi.iter);
fprintf('  relative error: %.1e\n', info_arnoldi.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_arnoldi.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_arnoldi.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_arnoldi.t_funm_quad);
fprintf('  thick restart time: %.1e\n', info_arnoldi.t_thick);

% sarnoldi.
rng(0);
param.thick = false;
param.arnoldi = 'sarnoldi';
[f_sarnoldi, info_sarnoldi] = funm_quad_log(A, b, param);
fprintf('sarnoldi\n');
fprintf('  time: %.1e\n', info_sarnoldi.t_total);
fprintf('  reason: %s\n', info_sarnoldi.stop_condition);
fprintf('  iteration number: %d\n', info_sarnoldi.iter);
fprintf('  relative error: %.1e\n', info_sarnoldi.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_sarnoldi.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_sarnoldi.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_sarnoldi.t_funm_quad);
fprintf('  thick restart time: %.1e\n', info_sarnoldi.t_thick);

% arnoldi (thick restart).
rng(0);
param.thick = true;
param.arnoldi = 'arnoldi';
[f_arnoldi_thick, info_arnoldi_thick] = funm_quad_log(A, b, param);
fprintf('arnoldi (thick restart)\n');
fprintf('  time: %.1e\n', info_arnoldi_thick.t_total);
fprintf('  reason: %s\n', info_arnoldi_thick.stop_condition);
fprintf('  iteration number: %d\n', info_arnoldi_thick.iter);
fprintf('  relative error: %.1e\n', info_arnoldi_thick.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_arnoldi_thick.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_arnoldi_thick.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_arnoldi_thick.t_funm_quad);
fprintf('  thick restart time: %.1e\n', info_arnoldi_thick.t_thick);

% sarnoldi (thick restart).
rng(0);
param.thick = true;
param.arnoldi = 'sarnoldi';
[f_sarnoldi_thick, info_sarnoldi_thick] = funm_quad_log(A, b, param);
fprintf('sarnoldi (thick restart)\n');
fprintf('  time: %.1e\n', info_sarnoldi_thick.t_total);
fprintf('  reason: %s\n', info_sarnoldi_thick.stop_condition);
fprintf('  iteration number: %d\n', info_sarnoldi_thick.iter);
fprintf('  relative error: %.1e\n', info_sarnoldi_thick.final_err);
fprintf('  matrix-vector products: %d\n', sum(info_sarnoldi_thick.num_oracle));
fprintf('  Arnoldi time: %.1e\n', info_sarnoldi_thick.t_arnoldi);
fprintf('  quadrature time: %.1e\n', info_sarnoldi_thick.t_funm_quad);
fprintf('  thick restart time: %.1e\n', info_sarnoldi_thick.t_thick);

% Error histories.
figure('Color', 'w');
semilogy(1:info_arnoldi.iter, info_arnoldi.err_history, '-o', ...
    'Color', [0, 0.4470, 0.7410], 'LineWidth', 2, ...
    'DisplayName', 'arnoldi');
hold on;
semilogy(1:info_sarnoldi.iter, info_sarnoldi.err_history, '-s', ...
    'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, ...
    'DisplayName', 'sarnoldi');
semilogy(1:info_arnoldi_thick.iter, info_arnoldi_thick.err_history, '--o', ...
    'Color', [0, 0.4470, 0.7410], 'LineWidth', 2, ...
    'DisplayName', 'arnoldi (thick restart)');
semilogy(1:info_sarnoldi_thick.iter, info_sarnoldi_thick.err_history, '--s', ...
    'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2, ...
    'DisplayName', 'sarnoldi (thick restart)');
legend('Location', 'best', 'Interpreter', 'none');
xlabel('Restart cycle');
ylabel('Relative error');
title('Logarithmic kernel: fixed with/without thick restart');

