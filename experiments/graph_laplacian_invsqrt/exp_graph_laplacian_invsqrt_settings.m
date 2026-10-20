clear;
close all;

matrix_data_path = "./data/";
if ~isfolder(matrix_data_path)
    mkdir(matrix_data_path);
end

if ispc || ismac
    data_path = "./data/";
elseif isunix
    data_path = "/scratch/jyliu/funm_quad_sketch_and_restart/frac_laplacian_invsqrt/data/";
end
if ~isfolder(data_path)
    mkdir(data_path);
end

figure_path = "./figure/";
if ~isfolder(figure_path)
    mkdir(figure_path);
end

if ispc || ismac
    graph_name = "p2p-Gnutella08";
else
    graph_name = "p2p-Gnutella25";
end

krylov_dim_list = [100, 150];

sfactor_list = 4;
fixed_method_list = ["arnoldi", "hmarnoldi", "sarnoldi", "shmarnoldi"];
adaptive_update_list = ["last_orth", "last_hmorth"];
thick_list = [false, true];
truncation_length_list = [2, 1, 0];

param = struct();
param.max_restarts = 500;
param.quad_tol = 1e-9;
param.max_num_quad_points = 1024;
param.stopping_accuracy = 1e-8;
param.number_thick = 5;
param.cond_tol = 1e6;


if ispc || ismac
    num_repeat = 1;
else
    num_repeat = 3;
end
seed_list = 0 : (num_repeat - 1);
