# Pull the container
docker pull nvcr.io/nvidia/cuda:12.4.1-devel-ubuntu22.04

# Enter the container and delete any modifications (--rm)
docker run --rm -it nvcr.io/nvidia/cuda:12.4.1-devel-ubuntu22.04 bash -l