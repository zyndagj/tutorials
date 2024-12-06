# Print out a message
#mpirun --allow-run-as-root echo hello
bcprun -no_redirect -c "echo hello"

# Print out the hostname of each process
#mpirun --allow-run-as-root hostname
bcprun -no_redirect -c "hostname"

# Run two processes per node
#mpirun -npernode 1 --allow-run-as-root hostname
bcprun -no_redirect -p 2 -c "hostname"

# Install a python package and load it
pip install lightning
bcprun -no_redirect -c "python test.py"

# Install the package on all nodes and then try loading it
bcprun -no_redirect -c "pip install lightning"
bcprun -no_redirect -c "python test.py"