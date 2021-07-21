#!/bin/bash
#SBATCH -o somd-array-gpu-%A.%a.out
#SBATCH -p GTX
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --time 24:00:00
#SBATCH --array=0-25

echo "CUDA DEVICES:" $CUDA_VISIBLE_DEVICES

lamvals=( 0.000 0.025 0.050 0.075 0.100 0.125 0.150 0.200 0.250 0.300 0.350 0.400 0.450 0.500 0.550 0.600 0.650 0.700 0.750 0.800 0.850 0.900 0.950 0.970 0.990 1.000 )
lam=${lamvals[SLURM_ARRAY_TASK_ID]}

if [ "$SLURM_ARRAY_TASK_ID" -eq "0" ]
then
    cat ../input/distres >> ../input/sim.cfg
fi

sleep 5

echo "lambda is: " $lam

mkdir lambda-$lam
cd lambda-$lam

export OPENMM_PLUGIN_DIR=/home/salome/sire.app/lib/plugins/

srun /home/salome/sire.app/bin/somd-freenrg -C ../../input/sim.cfg -l $lam -p CUDA
cd ..

wait

if [ "$SLURM_ARRAY_TASK_ID" -eq "25" ]
then
    sleep 30
    sbatch ../ljcor.sh
    sleep 30 
    sbatch ../standardstate.sh
    sleep 30 
    sbatch ../mbar.sh
fi

