# build the container
apptainer build nbody-efficient.sif Definition.nbody-efficient

# Look at image size
ls -lh nbody-efficient.sif