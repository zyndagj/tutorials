#!/bin/bash

# Time running training script on all nodes and GPUs in job
time bcprun -n ${NGC_ARRAY_SIZE} -p ${NGC_GPUS_PER_NODE} -no_redirect -c 'python pt_ddp_example.py'