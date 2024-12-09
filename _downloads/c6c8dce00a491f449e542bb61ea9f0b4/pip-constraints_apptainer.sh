# Delete old overlay and recreate
rm -rf ${APPTAINER_CACHEDIR}/pytorch_24.03_overlay
mkdir -p ${APPTAINER_CACHEDIR}/pytorch_24.03_overlay

# Launch a shell in the cuda devel container with
#   - fakeroot - appear to be root in the container
#   - overlay - allow modifications to be written to overlay directory
apptainer shell --fakeroot \
    --overlay ${APPTAINER_CACHEDIR}/pytorch_24.03_overlay \
    pytorch_24.03-py3.sif

# Save all existing packages and versions to a text file
pip list | awk '{print$1"=="$2}' | tail +3 > /root/base_constraints.txt

# Install any new packages without upgrading existing packages
pip install -c /root/base_constraints.txt torch torchvision torchaudio