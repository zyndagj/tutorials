# Your Docker hub username
HUB_USER=

# build the container
docker build -t ${HUB_USER}/nbody -f Dockerfile.nbody .

# Look at image size
docker images | grep nbody