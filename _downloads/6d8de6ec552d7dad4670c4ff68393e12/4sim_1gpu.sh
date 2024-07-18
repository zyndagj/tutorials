# Variable of shared arguments
GMX_ARGS="-ntmpi 1 -ntomp 8 -npme 0 -s benchMEM.tpr -cpt 1440 -nsteps 10000 -v -noconfout -nb gpu -dlb yes -bonded gpu"

# without MPS
# unique log and ouptput files
# processes are backgrounded and stdout sent to /dev/null
gmx mdrun ${GMX_ARGS} -gpu_id 0 -e /tmp/out1 -g 1gpu_sim1-4.log &> /dev/null &
gmx mdrun ${GMX_ARGS} -gpu_id 0 -e /tmp/out2 -g 1gpu_sim2-4.log &> /dev/null &
gmx mdrun ${GMX_ARGS} -gpu_id 0 -e /tmp/out3 -g 1gpu_sim3-4.log &> /dev/null &
gmx mdrun ${GMX_ARGS} -gpu_id 0 -e /tmp/out4 -g 1gpu_sim4-4.log &> /dev/null &
# wait for processes to complete
wait

# with MPS
nvidia-cuda-mps-control -d
# unique log and ouptput files
# processes are backgrounded and stdout sent to /dev/null
gmx mdrun ${GMX_ARGS} -gpu_id 0 -e /tmp/out1 -g 1gpu_sim1-4_mps.log &> /dev/null &
gmx mdrun ${GMX_ARGS} -gpu_id 0 -e /tmp/out2 -g 1gpu_sim2-4_mps.log &> /dev/null &
gmx mdrun ${GMX_ARGS} -gpu_id 0 -e /tmp/out3 -g 1gpu_sim3-4_mps.log &> /dev/null &
gmx mdrun ${GMX_ARGS} -gpu_id 0 -e /tmp/out4 -g 1gpu_sim4-4_mps.log &> /dev/null &
# wait for processes to complete
wait
# stop MPS server
echo quit | nvidia-cuda-mps-control
