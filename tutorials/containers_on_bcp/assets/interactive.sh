#!/bin/bash

# Create a workspace
ngc workspace create --name ${USER}_containers

# Launch job and mount workspace
ngc batch run --name "python-devel" --result /results \
	--total-runtime 1h --instance dgxa100.80g.1.norm \
        --commandline "jupyter lab --ip=0.0.0.0 --allow-root --no-browser --NotebookApp.token='' --notebook-dir=/ --NotebookApp.allow_origin='*'" \
        --port 8888 --workspace ${USER}_containers:/mnt/workspace:RW \
	--image "nvcr.io/nvidia/pytorch:24.03-py3"
