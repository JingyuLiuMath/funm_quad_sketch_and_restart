clear;
close all;

data_path = "./data/";
figure_path = "./figure/";
if ~isfolder(data_path)
    mkdir(data_path);
end
if ~isfolder(figure_path)
    mkdir(figure_path);
end

% Use the full QCD example on every platform.
sfactor_list = [2, 4];
fixed_method_list = ["arnoldi", "hmarnoldi", "sarnoldi", "shmarnoldi"];
adaptive_update_list = ["last_orth", "last_hmorth", "last_sorth", "last_shmorth"];
thick_list = [false, true];
if ispc
    num_repeat = 1;
else
    num_repeat = 3;
end
seed_list = 0 : (num_repeat - 1);
truncation_length_list = [2, 1, 0];

param = struct();
if ispc
    param.krylov_dim = 150;
else
    param.krylov_dim = 150;
end
param.max_restarts = 300;
param.quad_tol = 1e-11;
param.max_num_quad_points = 1024;
param.stopping_accuracy = 1e-10;
param.number_thick = 5;
param.cond_tol = 1e6;
