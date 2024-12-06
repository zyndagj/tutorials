# Create an overlay directory
#   - The base container is never changed, just the overlay
#   - Overlay only works with a single image
mkdir -p ${APPTAINER_CACHEDIR}/cuda-devel_overlay

# Launch a shell in the cuda devel container with
#   - fakeroot - appear to be root in the container
#   - overlay - allow modifications to be written to overlay directory
apptainer shell --fakeroot \
    --overlay ${APPTAINER_CACHEDIR}/cuda-devel_overlay \
    cuda_12.4.1-devel-ubuntu22.04.sif

# Your cluster may also support overlay images instead of directories