# You should already be on a compute node
# srun -p gpu -n 1 -G 1 --cpus-per-task 16 -t 03:00:00 --pty bash -l

# build the container
apptainer build nbody.sif Definition.nbody

# Look at image size
ls -lh nbody.sif