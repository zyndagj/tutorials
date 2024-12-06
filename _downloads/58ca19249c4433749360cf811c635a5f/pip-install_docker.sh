# Launch pytorch container on your host
docker run --rm -it nvcr.io/nvidia/pytorch:24.03-py3 bash -l

# Pip install pytorch inside the running container
pip install torch torchvision torchaudio