# Run nbody benchmark
#   - Include GPU support (--nv)
#   - Exclude $HOME mount (-c)
#   - Mount CWD (-B)
apptainer shell --nv -c -B $PWD:$PWD lightning.sif