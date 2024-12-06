#!/bin/bash
#SBATCH --job-name=ptl_ddp_example
#SBATCH --nodes=2             # Set number of nodes
#SBATCH --tasks-per-node=2    # Set number of GPUs per node
#SBATCH --gpus-per-node=2     #   - set to the same
#SBATCH --mem=16GB            # Set memory limits (consider --exclusive)
#SBATCH --cpus-per-task=8
#SBATCH --output=%x-%j.out
#SBATCH --partition=gpu
#SBATCH --time=00:30:00

#export LOGLEVEL=INFO
#export NCCL_DEBUG=INFO

export MASTER_ADDR=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n 1)
export MASTER_PORT=$(expr 10000 + $(echo -n ${SLURM_JOBID} | tail -c 4))

# Launches one container per node
#  - Contaienr spawns multiple processes
srun --mpi=pmi2 apptainer exec --nv lightning.sif \
	bash -c "export NODE_RANK=\${SLURM_PROCID}; python ptl_ddp_example.py -N $SLURM_JOB_NUM_NODES -p ${SLURM_GPUS_ON_NODE}"
