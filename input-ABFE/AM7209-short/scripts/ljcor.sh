#!/bin/bash

cd output/lambda-0.000
~/sire.app/bin/lj-tailcorrection -C ../../input/sim.cfg -l 0.00 -r traj000000001.dcd -s 1 1> ../freenrg-LJCOR-lam-0.000.dat 2> /dev/null

cd ../..
cd output/lambda-1.000
~/sire.app/bin/lj-tailcorrection -C ../../input/sim.cfg -l 1.00 -r traj000000001.dcd -s 1 1> ../freenrg-LJCOR-lam-1.000.dat 2> /dev/null
cd ../..

wait

# utility script to get final LJ correction term
python parselj.py freenrg-LJCOR-lam-0.000.dat freenrg-LJCOR-lam-1.000.dat > freenrg-LJCOR.dat
