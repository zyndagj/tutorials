# Create an overlay directory
#   - The base container is never changed, just the overlay
#   - Overlay only works with a single image
mkdir -p ${APPTAINER_CACHEDIR}/pytorch_24.03_overlay

# Launch a shell in the cuda devel container with
#   - fakeroot - appear to be root in the container
#   - overlay - allow modifications to be written to overlay directory
apptainer shell --fakeroot \
    --overlay ${APPTAINER_CACHEDIR}/pytorch_24.03_overlay \
    pytorch_24.03-py3.sif

# Pip install pytorch inside the running container
pip install torch torchvision torchaudio