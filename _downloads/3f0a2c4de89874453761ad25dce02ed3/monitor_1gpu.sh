# Kill child processes on exit
trap 'pkill -P $$' SIGINT SIGTERM EXIT

# Start capturing utilization to CSV file in the background
nvidia-smi -i 0 --query-gpu=timestamp,name,utilization.gpu,memory.used --format=csv -l 5 > utilization.csv &

# Runs benchMEM with 1 process, 8 threads, for 10k steps on GPU 0
# Simulation output is written to /tmp/out (deleted)
# Log is written to 1gpu.log
gmx mdrun -ntmpi 1 -ntomp 8 -npme 0 -s benchMEM.tpr -cpt 1440 \
	-nsteps 10000 -v -noconfout -nb gpu -dlb yes -gpu_id 0 \
	-bonded gpu -e /tmp/out -g monitor_1gpu.log
