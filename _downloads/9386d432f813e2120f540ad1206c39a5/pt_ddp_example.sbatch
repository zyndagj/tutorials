#!/bin/bash
#SBATCH --job-name=pt_ddp_example
#SBATCH --nodes=2             # Set number of nodes
#SBATCH --gpus-per-node=2     # Set number of GPUs per node
#SBATCH --mem=32GB            # Set memory limits (consider --exclusive)
#SBATCH --tasks-per-node=1
#SBATCH --output=%x-%j.out
#SBATCH --cpus-per-task=8
#SBATCH --partition=gpu
#SBATCH --time=00:30:00

# Job debug info
echo "Launching on ${SLURM_JOB_NUM_NODES} nodes"
echo "Launching on: " ${SLURM_JOB_NODELIST}
echo "Launching ${SLURM_NTASKS_PER_NODE} tasks per node"
echo "Using ${SLURM_GPUS_ON_NODE} GPUs per task"

# Optional debug logging
#export LOGLEVEL=INFO
#export NCCL_DEBUG=INFO

##### No need to edit these #########################################################
# main address is detected by first name in nodelist
export MASTER_ADDR=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n 1)
# port is chosen by jobID (prevents collisions if nodes are shared)
export MASTER_PORT=$(expr 10000 + $(echo -n ${SLURM_JOBID} | tail -c 4))
export WORLD_SIZE=$((${SLURM_NNODES} * ${SLURM_GPUS_ON_NODE}))
#####################################################################################
echo "Training on ${WORLD_SIZE} GPUs - ${MASTER_ADDR}:${MASTER_PORT}"

srun --mpi=pmi2 apptainer exec --nv lightning.sif torchrun \
	--nnodes ${SLURM_JOB_NUM_NODES} \
	--nproc_per_node ${SLURM_GPUS_ON_NODE} \
	--rdzv_id $RANDOM \
	--rdzv_backend c10d \
	--rdzv_endpoint $MASTER_ADDR:${MASTER_PORT} \
	pt_ddp_example.py
