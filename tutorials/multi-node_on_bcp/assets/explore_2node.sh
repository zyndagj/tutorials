# Print out a message
mpirun --allow-run-as-root echo hello

# Print out the hostname of each process
mpirun --allow-run-as-root hostname

# Run two processes per node
mpirun -npernode 2 --allow-run-as-root hostname

# Install a python package and load it
pip install lightning
mpirun -npernode 2 --allow-run-as-root python -c 'import pytorch_lightning'
mpirun -npernode 2 --allow-run-as-root python -c 'import pytorch_lightning'
bcprun -no_redirect -c "pip install biopython"