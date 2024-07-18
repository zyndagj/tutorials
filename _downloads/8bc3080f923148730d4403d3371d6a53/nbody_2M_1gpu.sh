#!/bin/bash

ngc batch run --name "nbody-test-2M-1gpu" \
	--total-runtime 1200s --instance dgxa100.80g.1.norm \
	--commandline "nbody -benchmark -numbodies=2000000" \
	--result /results \
	--image "nvcr.io/${ORG}${TEAM}/${USER}_nbody:efficient"
