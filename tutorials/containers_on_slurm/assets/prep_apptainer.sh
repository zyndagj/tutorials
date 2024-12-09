# Start a 3 hour interactive job with 1 GPU
srun -p interactive -n 1 -G 1 --cpus-per-task 16 -t 03:00:00 --pty bash -l

# Change cachedir to /tmp
export APPTAINER_CACHEDIR=/tmp/${USER}_apptainer_cache

# Create workspace for tutorial
mkdir -p ${MYDATA}/containers

# Change to your workspace
#  - This is a good place to store containers
#  - Keep definition files on $HOME
cd ${MYDATA}/containers

# Pull the container
apptainer pull docker://nvcr.io/nvidia/pytorch:24.03-py3