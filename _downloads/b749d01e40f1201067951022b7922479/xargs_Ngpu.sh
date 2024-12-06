#!/bin/bash

# Number of cores in job
export NCORES=22
# Number of tasks to run concurrently
export NP=${1:-8}
# Number of tasks in queue
export NT=16
# Number of GPUs
export NG=2

# Calculate total number of concurrent tasks
export TT=$(( $NP * $NG ))
# Variable of shared arguments
export GMX_ARGS="-ntmpi 1 -ntomp $(( $NCORES / $NP )) -npme 0 -s benchMEM.tpr -cpt 1440 -nsteps 10000 -v -noconfout -nb gpu -dlb yes -bonded gpu"

# Define and export function for running GMX
function run_gmx {
	TID=$1
	GPU=$(( $SLOT % $NG ))
	# Load GMX environment
	. /usr/local/gromacs/avx2_256/bin/GMXRC.bash
	echo starting task ${TID} in slot ${SLOT} on GPU ${GPU}
	# Task is not backgrounded
	gmx mdrun ${GMX_ARGS} -gpu_id ${GPU} -e /tmp/out${TID} -g ${NG}gpu_sim${TID}-${TT}_xargs.log &> /dev/null
	echo finished task ${TID} in slot ${SLOT}
}
export -f run_gmx

# with MPS
echo "Running ${NP} simulations concurrently with MPS"
nvidia-cuda-mps-control -d
# Spawn NP processes and wait for them to complete
for i in $(seq 1 $NT); do
	echo run_gmx $i
done | xargs -P ${TT} --process-slot-var=SLOT -I {} bash -c '{}'
# stop MPS server
echo quit | nvidia-cuda-mps-control
