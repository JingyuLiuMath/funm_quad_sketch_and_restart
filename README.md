# Sketch-and-Restart for Matrix Functions

MATLAB implementations and scripts for reproducing the numerical experiments associated with:

> Stefan Güttel, Jingyu Liu, and Lauri Nyman,
> **Sketch-and-Restart: Randomized Sketching in Quadrature-Based Restarting for Matrix Functions**,
> [arXiv:2607.10354](https://arxiv.org/abs/2607.10354), 2026.

The repository combines quadrature-based restarting with Arnoldi and harmonic Arnoldi methods, their sketched variants, and adaptive methods based on truncated Arnoldi. It includes comparisons of convergence, runtime, storage, quadrature costs, and thick restarting. The current experiment configurations have evolved since arXiv v1; use the settings files and matching input data when reproducing a particular set of results.

## Repository structure

```text
src/
  funm_quad/                    Matrix-function solvers and quadrature routines
  arnoldi/                      Standard, harmonic, sketched, and adaptive Arnoldi
  thick_restarting/             Schur-based subspace retention
  tools/                        Sparse sign sketches
  matrices/                     Laplacian and convection-diffusion utilities
test/
  invsqrt/                      Inverse-square-root examples
  log/                          Logarithmic-function examples
experiments/
  qcd_invsqrt/                   QCD inverse-square-root experiments
  graph_laplacian_invsqrt/       Fractional graph Laplacian experiments
  conv_diff_log/                 Convection-diffusion logarithm experiments
sbatch/                         Slurm submission scripts for the experiments
setuppath.m                     Add source, tests, and experiments to the path
```

## Requirements and setup

Use MATLAB with support for string arrays and `exportgraphics`. The supplied Slurm scripts request MATLAB R2023b; Slurm and environment modules are needed only for cluster submission.

Open MATLAB in the repository root and run:

```matlab
setuppath;
```

`setuppath` resolves the repository location from its own file. Run each experiment from its corresponding directory so the relative input and output paths resolve correctly.

## Quick start

From the repository root:

```matlab
setuppath;
test_invsqrt_fix;       % Compare four fixed-dimension methods
% test_invsqrt_adaptive;
% test_log_fix;
% test_log_adaptive;
```

These examples construct a two-dimensional Dirichlet Laplacian, compute a spectral reference solution, and display relative-error histories. The `*_thick_restart` scripts in the same test directories compare ordinary and thick restarting.

## Solver inputs and outputs

The two entry points are:

```matlab
[f, info] = funm_quad_invsqrt(A, b, param);
[g, info] = funm_quad_log(A, b, param);
```

`funm_quad_invsqrt` approximates `A^(-1/2)*b`. `funm_quad_log` approximates `g(A)*b`, where `g(z) = log(1+z)/z`; for nonsingular `A`, this is `log(I+A)*(A\b)`. The inverse square root uses the principal branch. The matrix and quadrature shifts must satisfy the assumptions of the corresponding method in the paper.

`A` can be a square matrix or a function handle computing `A*v`, and `b` is a column vector. Set the parameters explicitly; the solvers do not supply defaults. Complete examples are provided in `test/invsqrt/` and `test/log/`.

| Parameter | Meaning |
| --- | --- |
| `mode` | `"fixed"` or `"adaptive"` |
| `arnoldi` | Fixed method: `"arnoldi"`, `"hmarnoldi"`, `"sarnoldi"`, or `"shmarnoldi"`; use `"arnoldi"` for adaptive methods |
| `krylov_dim`, `max_restarts` | Maximum new Arnoldi steps per cycle and maximum number of cycles; adaptive runs reduce the step limit to the first-cycle dimension |
| `truncation_length` | Number of preceding basis vectors used in truncated orthogonalization; fixed experiments use `inf` |
| `sketch_dim` | Sketch dimension for fixed sketched methods |
| `sketching_size_control`, `cond_tol` | Sketch-size control and basis-condition threshold for adaptive methods |
| `update` | Adaptive last-vector correction: `"last_orth"`, `"last_hmorth"`, `"last_sorth"`, or `"last_shmorth"` |
| `thick`, `number_thick` | Enable thick restarting and choose the requested number of retained modes |
| `quad_tol`, `max_num_quad_points` | Quadrature tolerance and node limit |
| `stopping_accuracy` | Relative correction threshold |
| `exact` | Optional reference vector for recording relative errors |

Fixed harmonic methods (`hmarnoldi`, `shmarnoldi`) use `thick=false`. Thick restarting is available for fixed Arnoldi/sketched Arnoldi and adaptive methods. Retained modes add to the new Arnoldi steps in a cycle.

The stopping test compares the correction norm with the norm of the updated approximation. Supplying `param.exact` records errors without changing this test.

| `info` field | Meaning |
| --- | --- |
| `stop_condition`, `iter` | Reason for stopping and number of cycles performed |
| `err_history`, `final_err` | Relative errors against `param.exact`; `[]` and `NaN` when no reference is supplied |
| `dim`, `num_thick` | Projected dimensions and retained-mode counts by cycle |
| `num_oracle` | Operator applications by cycle; sum for the total count |
| `num_quadpoints`, `quad_converged` | Quadrature node counts and convergence flags |
| `update_norm` | Absolute correction norms |
| `t_total`, other `t_*`, `arnoldi_info` | Total runtime and operation-level timings |

For QCD, one operator application evaluates `Q*(Q*v)` and therefore involves two multiplications by `Q`.

## Reproducing the numerical experiments

The three experiment directories correspond to the following problems in the paper:

| Directory | Matrix-function action |
| --- | --- |
| `experiments/qcd_invsqrt/` | `(Q^2)^(-1/2)*b`, with `b = Q*e1/norm(Q*e1)`, giving the normalized matrix-sign action |
| `experiments/graph_laplacian_invsqrt/` | `A^(-1/2)*b`, with graph Laplacian `A` and `b = A*b0` |
| `experiments/conv_diff_log/` | `g(A)*b`, where `g(z) = log(1+z)/z`, for a two-dimensional convection-diffusion matrix |

Edit the corresponding `exp_*_settings.m` before running. These files select dimensions, tolerances, methods, sketch factors, truncation lengths, repetitions, and paths. They clear workspace variables, so configure them in the files. Each experiment runs fixed and adaptive methods with and without thick restarting where supported.

Current settings use sketch factor `4`, adaptive truncation lengths `[2, 1, 0]`, and five retained modes. QCD and graph experiments use fixed dimensions `[100, 150]`; convection-diffusion uses `[150, 200]`. Windows/macOS use one random seed; Linux uses three. Plot and table scripts also contain selections such as `selection.krylov_dim`; keep these consistent with the runs being reproduced.

### Input data and output paths

The experiment scripts load matrices, right-hand sides, and reference vectors from their own `./data/` directories:

| Experiment | Input files |
| --- | --- |
| QCD | `qcd_matrix.mat` (`Q`, `b`) and `qcd_sol.mat` (`f_ex`) |
| Graph Laplacian | `frac_laplacian_matrix_<graph_name>.mat` (`A`, `b`) and `frac_laplacian_sol_<graph_name>.mat` (`f_ex`) |
| Convection-diffusion | `conv_diff_n_<n>_nu_<nu>_matrix.mat` (`A`, `b`) and the matching `_sol.mat` (`f_ex`) |

The supplied graph reference files cover `p2p-Gnutella08` and `p2p-Gnutella31`; graph settings select the former on Windows/macOS and the latter on Linux. The supplied convection-diffusion data use `n=512`, `nu=0.1`; set `n=512` to use these files on Windows/macOS, whose default is `n=256`.

The graph and convection-diffusion `exp_*_ref.m` files contain reference-generation code that is currently commented out. Enable and configure that code if new input/reference files are needed; running the files as supplied does not generate missing data.

On Windows/macOS, results are saved in `./data/`. Linux settings use paths under `/scratch/jyliu/funm_quad_sketch_and_restart/`; change `data_path` for your machine. Figures are exported as PDF and PNG to `./figure/`.

### QCD

From the repository root after `setuppath`:

```matlab
cd('experiments/qcd_invsqrt');
exp_qcd_invsqrt_fix;
exp_qcd_invsqrt_adaptive;
print_qcd_invsqrt;
plot_qcd_invsqrt_err;
plot_qcd_invsqrt_total_ratios;
plot_qcd_invsqrt_arnoldi_ratios;
cd('../..');
```

### Fractional graph Laplacian

```matlab
cd('experiments/graph_laplacian_invsqrt');
exp_graph_laplacian_invsqrt_fix;
exp_graph_laplacian_invsqrt_adaptive;
print_graph_laplacian_invsqrt;
plot_graph_laplacian_invsqrt_err;
plot_graph_laplacian_invsqrt_total_ratios;
plot_graph_laplacian_invsqrt_arnoldi_ratios;
print_graph_laplacian_invsqrt_thick_restarting;
plot_graph_laplacian_invsqrt_thick_restarting;
plot_graph_laplacian_invsqrt_quadpoints;
cd('../..');
```

### Convection-diffusion

```matlab
cd('experiments/conv_diff_log');
exp_conv_diff_log_fix;
exp_conv_diff_log_adaptive;
print_conv_diff_log;
plot_conv_diff_log_err;
plot_conv_diff_log_total_ratios;
plot_conv_diff_log_arnoldi_ratios;
print_conv_diff_log_thick_restarting;
plot_conv_diff_log_thick_restarting;
cd('../..');
```

The experiments save per-run approximations, errors, parameters, seeds, and timings in MAT files. Postprocessing selects the run with median total runtime for each configuration. The error plots use reference-based relative errors, and the timing plots separate total solver costs from Arnoldi costs.

For cluster runs, adapt the repository paths, node requests, and MATLAB module in `sbatch/`, then submit the corresponding fixed/adaptive job scripts with `sbatch`.

## License

The source code is distributed under the [MIT License](LICENSE).
