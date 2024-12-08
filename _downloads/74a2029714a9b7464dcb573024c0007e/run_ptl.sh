#!/bin/bash

# Set address for main process
export NGC_MASTER_ADDR=launcher-svc-${NGC_JOB_ID}

# Launch PTL script on all nodes and GPUs in job
bcprun -no_redirect -n ${NGC_ARRAY_SIZE} -p ${NGC_GPUS_PER_NODE} -c "python ptl_ddp_example.py -N ${NGC_ARRAY_SIZE} -p ${NGC_GPUS_PER_NODE}"