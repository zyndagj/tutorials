#!/bin/bash

# Number of cores in job
NCORES=22
# Number of tasks to run concurrently
NP=${1:-8}

# Variable of shared arguments
GMX_ARGS="-ntmpi 1 -ntomp $(( $NCORES / $NP )) -npme 0 -s benchMEM.tpr -cpt 1440 -nsteps 10000 -v -noconfout -nb gpu -dlb yes -bonded gpu"

# without MPS
echo "Running ${NP} simulations concurrently without MPS"
# Spawn NP processes and wait for them to complete
for i in $(seq 1 $NP); do
	gmx mdrun ${GMX_ARGS} -gpu_id 0 -e /tmp/out${i} -g 1gpu_sim${i}-${NP}.log &> /dev/null &
done
# wait for processes to complete
wait

# with MPS
echo "Running ${NP} simulations concurrently with MPS"
nvidia-cuda-mps-control -d
# Spawn NP processes and wait for them to complete
for i in $(seq 1 $NP); do
	gmx mdrun ${GMX_ARGS} -gpu_id 0 -e /tmp/out${i} -g 1gpu_sim${i}-${NP}_mps.log &> /dev/null &
done
# wait for processes to complete
wait
# stop MPS server
echo quit | nvidia-cuda-mps-control
