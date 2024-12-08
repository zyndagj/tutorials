# Run on 2 GPUs of any type
#  (-g) argument sets how many GPUs each process will use
srun -p gpu -N 2 -n 2 --gpus-per-node 1 --mpi=pmi2 apptainer exec --nv lightning.sif all_reduce_perf_mpi -b 1G -e 4G -f 2 -g 1

# Run on 4 H100 GPUs across 2 nodes
srun -p gpu --mem=32G -N 2 -n 2 --gpus-per-node h100:2 --mpi=pmi2 apptainer exec --nv lightning.sif all_reduce_perf_mpi -b 1G -e 4G -f 2 -g 2