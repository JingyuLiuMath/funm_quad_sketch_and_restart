#!/bin/bash

#SBATCH --job-name=conv_diff_log_fix
#SBATCH --output=conv_diff_log_fix_%j.out
#SBATCH --error=conv_diff_log_fix_%j.err
#SBATCH --nodelist=bigMem1
#SBATCH --exclusive
#SBATCH --time=18:00:00

module unload MATLAB
module load MATLAB/R2023b || exit $?

echo "=========================================="
echo "        SLURM JOB INFORMATION"
echo "=========================================="
echo "Job ID:           $SLURM_JOB_ID"
echo "Job Name:         $SLURM_JOB_NAME"
echo "Submit Directory: $SLURM_SUBMIT_DIR"
echo "Submit Host:      $SLURM_SUBMIT_HOST"
echo "Partition:        $SLURM_JOB_PARTITION"
echo "Node List:        $SLURM_JOB_NODELIST"
echo "Nodes:            $SLURM_NNODES"
echo "Tasks:            $SLURM_NTASKS"
echo "CPUs per Task:    $SLURM_CPUS_PER_TASK"
echo "Total CPUs:       $SLURM_NPROCS"
echo "Working Directory: $(pwd)"
echo "=========================================="
echo ""
echo "        RESOURCE ALLOCATION"
echo "=========================================="
echo "Memory per Node:  $SLURM_MEM_PER_NODE"
echo "Time Limit:       $SLURM_TIMELIMIT"
echo "=========================================="
echo ""
echo "        ENVIRONMENT VARIABLES"
echo "=========================================="
echo "PATH:             $PATH"
echo "LD_LIBRARY_PATH:  $LD_LIBRARY_PATH"
echo "LIBRARY_PATH:     $LIBRARY_PATH"
echo "CPATH:            $CPATH"
echo "=========================================="
echo ""
echo "        MODULE INFORMATION"
echo "=========================================="
module list 2>&1
echo "=========================================="
echo ""
echo "Job started at: $(date)"
echo "=========================================="
echo ""

matlab -batch "cd('/home/jyliu/funm_quad_sketch_and_restart'); setuppath; cd('experiments/conv_diff_log'); exp_conv_diff_log_fix;"
matlab_status=$?

echo ""
echo "MATLAB finished at: $(date)"
echo "MATLAB exit status: $matlab_status"
echo "=========================================="
if [ "$matlab_status" -eq 0 ]; then
    echo "Job completed"
else
    echo "Job failed"
fi
echo "=========================================="
exit "$matlab_status"
