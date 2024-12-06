# Your Docker hub username
HUB_USER=

# build the container
docker build -t ${HUB_USER}/lighting:latest -f Dockerfile.lightning .

# Push the container
docker push ${HUB_USER}/nbody:efficient