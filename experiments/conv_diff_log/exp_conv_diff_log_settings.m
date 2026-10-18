clear;
close all;

if ispc || ismac
    data_path = "./data/";
elseif isunix
    data_path = "/scratch/jyliu/funm_quad_sketch_and_restart/conv_diff_log/data/";
end
figure_path = "./figure/";
if ~isfolder(data_path)
    mkdir(data_path);
end
if ~isfolder(figure_path)
    mkdir(figure_path);
end

figure_prefix = "./figure/conv_diff_log";

if ispc
    n = 50;
else
    n = 180;
end
N = n^2;
nu = 1;
sfactor_list = [2, 4];
fixed_method_list = ["arnoldi", "hmarnoldi", "sarnoldi", "shmarnoldi"];
adaptive_update_list = ["last_orth", "last_hmorth", "last_sorth", "last_shmorth"];
thick_list = [false, true];
num_repeat = 3;
seed_list = 0 : (num_repeat - 1);
truncation_length_list = [2, 1, 0];

param = struct();
if ispc
    param.krylov_dim = 30;
else
    param.krylov_dim = 150;
end
param.max_restarts = 30;
param.quad_tol = 1e-11;
param.max_num_quad_points = 512;
param.stopping_accuracy = 1e-10;
param.number_thick = 5;
param.cond_tol = 1e6;
