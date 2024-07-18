# Runs benchMEM with 2 processes, 8 threads each (16 total), for 10k steps on GPUs 0 and 1
# Simulation output is written to /dev/null (deleted)
# Log is written to 2gpu.log
gmx mdrun -ntmpi 2 -ntomp 8 -npme 0 -s benchMEM.tpr -cpt 1440 \
	-nsteps 10000 -v -noconfout -nb gpu -dlb yes -gpu_id 01 \
	-bonded gpu -e /tmp/out -g 2gpu.log
