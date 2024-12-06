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

##### Number of total processes 
echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "
echo "Nodelist:= " ${SLURM_JOB_NODELIST}
echo "Number of nodes:= " ${SLURM_JOB_NUM_NODES}
echo "Ntasks per node:= " ${SLURM_GPUS_ON_NODE}
echo "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX "

# ******************* These are read internally it seems ***********************************
# ******** Master port, address and world size MUST be passed as variables for DDP to work 
#export LOGLEVEL=INFO
#export NCCL_DEBUG=INFO

export MASTER_ADDR=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | head -n 1)
export MASTER_PORT=$(expr 10000 + $(echo -n ${SLURM_JOBID} | tail -c 4))
export WORLD_SIZE=$((${SLURM_NNODES} * ${SLURM_GPUS_ON_NODE}))
echo "Training on ${WORLD_SIZE} GPUs - ${MASTER_ADDR}:${MASTER_PORT}"
# ******************************************************************************************

srun --mpi=pmi2 apptainer exec --nv lightning.sif torchrun \
	--nnodes ${SLURM_JOB_NUM_NODES} \
	--nproc_per_node ${SLURM_GPUS_ON_NODE} \
	--rdzv_id $RANDOM \
	--rdzv_backend c10d \
	--rdzv_endpoint $MASTER_ADDR:${MASTER_PORT} \
	pt_ddp_example.py
