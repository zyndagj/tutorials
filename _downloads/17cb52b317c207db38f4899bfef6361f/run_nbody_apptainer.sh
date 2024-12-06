# Check that GPU is still detected
apptainer exec --nv nbody-efficient.sif nvidia-smi

# Run nbody benchmark
apptainer exec --nv nbody-efficient.sif nbody -benchmark -numbodies=2000000