mkdir 
[ -e out ] && rm -rf out                                                                 
mkdir out && cd out                                                                      

# one gpu
gmx mdrun -ntmpi 1 -ntomp 8 -npme 0 -s /raid/benchMEM.tpr -cpt 1440 -nsteps 10000 -v -noconfout -nb gpu -dlb yes -gpu_id 0 -bonded gpu -e /dev/null -g 1gpu.log

# two gpus
gmx mdrun -ntmpi 2 -ntomp 8 -npme 0 -s /raid/benchMEM.tpr -cpt 1440 -nsteps 10000 -v -noconfout -nb gpu -dlb yes -gpu_id 01 -bonded gpu -e /dev/null -g 2gpu.log
