# Start a 2 hour interactive job with 1 GPU
srun -p gpu -n 1 -G 1 --cpus-per-task 16 -t 03:00:00 --pty bash -l

# Pull the container
enroot import docker://nvcr.io#nvidia/cuda:12.4.1-devel-ubuntu22.04