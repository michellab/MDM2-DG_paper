#!/bin/bash
#SBATCH -o analyse-free-nrg-%A.%a.out
#SBATCH -p GTX
#SBATCH -n 1
#SBATCH --time 01:00:00

#sleep 30 

export OPENMM_PLUGIN_DIR=/home/salome/sire.app/lib/plugins/

#srun /home/$USER/sire.app/bin/analyse_freenrg mbar -i 'lambda*/simfile*' -p 95 --overlap --temperature 298.0 > freenrg-MBAR-p-95-overlap.dat
srun /home/salome/sire.app/bin/analyse_freenrg mbar -i output/lambda*/simfile.dat -p 95 --overlap --temperature 298.0 > freenrg-MBAR-p-95-overlap.dat

bzip2 lambda-*/*




