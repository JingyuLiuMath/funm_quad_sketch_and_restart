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

if ispc
    n = 50;
else
    n = 150;
end
D = 1e-2;
sfactor_list = [4, 8];
fixed_method_list = ["arnoldi", "hmarnoldi", "sarnoldi", "shmarnoldi"];
adaptive_update_list = ["last_orth", "last_hmorth", "last_sorth", "last_shmorth"];
thick_list = [false, true];
num_repeat = 3;
seed_list = 0 : (num_repeat - 1);
truncation_length_list = [2, 1, 0];

param = struct();
param.krylov_dim = 30;
param.max_restarts = 300;
param.quad_tol = 1e-9;
param.max_num_quad_points = 512;
param.stopping_accuracy = 1e-8;
param.number_thick = 5;
param.sketching_mat_type = "sparse sign";
param.cond_tol = 1e6;
