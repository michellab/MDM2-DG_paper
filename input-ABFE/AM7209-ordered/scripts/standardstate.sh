#!/bin/bash
#SBATCH -o standardstate-%A.%a.out
#SBATCH -p GTX
#SBATCH -n 1
#SBATCH --time 03:00:00

sleep 30 

cd lambda-1.000

srun /home/$USER/sire.app/bin/standardstatecorrection -C ../../input/sim.cfg -r traj000000001.dcd -s 1 -b 4 -d 0.25 -o 6 > ../standard-state-s-1-b-4-d-0.25-o-6.dat

cd ..

wait


