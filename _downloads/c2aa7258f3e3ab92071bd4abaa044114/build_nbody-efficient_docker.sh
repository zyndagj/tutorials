# Your Docker hub username
HUB_USER=

# build the container
docker build -t ${HUB_USER}/nbody:efficient -f Dockerfile.nbody-efficient .

# Look at image size
docker images | grep nbody

# Push the container
docker push ${HUB_USER}/nbody:efficient