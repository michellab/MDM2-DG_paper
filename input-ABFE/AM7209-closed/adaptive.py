#!/usr/bin/env python
# coding: utf-8

# 1) Work out binding free energies using all available sampling in the subfolders
# 2) Work out list of windows whose stdev exceeds threshold 
# 3) Generate submission scripts to carry out one round of adaptive sampling
import glob,sys,os
import copy
import numpy as np
import math
import matplotlib.pyplot as plt
### ADAPTIVE SAMPLING ROUTE
NRUNS = 5
#NEPOCHS = 12
EPOCHTIME = 5 # ns  
DISCARD = 0.5 # ns
CHUNKSIZE = 0.5 # ns
THRESHOLD = 0.1 #kcal/mol
TIMESTEP = 0.000004 #ns
NRGFREQ = 250
### HELPER FUNCTIONS
def doMBAR(simfiles, startime, endtime):
    cmd = "rm -rf tmp ; mkdir tmp" 
    os.system(cmd)
    cumtime = 0.0
    for simfile in simfiles:
        #print (simfile, max_time)
        # Ugly
        lamval = os.path.split(simfile)[-2].split("/")[-1]
        #print (lamval)
        istream = open(simfile,'r')
        ostream = open("tmp/%s.dat" % lamval,'w')
        ifile = istream.readlines()
        for line in ifile:
            if line.startswith("#"):
                ostream.write(line)
                continue
            elems = line.split()
            time = float(elems[0])*TIMESTEP
            if (time < startime):
                continue
            if (time > endtime):
                break
            ostream.write(line)
        cumtime += (time - startime)
        istream.close()
        ostream.close()
    # TODO ADD OVERLAP MATRIX TO OUTPUT
    cmd = "~/sire.app/bin/analyse_freenrg mbar -i tmp/*.dat -o tmp/mbar.dat -p 100"
    os.system(cmd)
    istream = open("tmp/mbar.dat","r")
    results = {}
    ofile = istream.readlines()
    inDGsection = False
    for line in ofile:
        if line.startswith("#DG from neighbouring"):
            inDGsection = True
            continue
        if line.startswith("#PMF from MBAR"):
            inDGsection = False
            continue
        if inDGsection:
            elems = line.split()
            lami = float(elems[0])
            lamj = float(elems[1])
            DGij = float(elems[2])
            DGij_sig = float(elems[3])
            results["%.5f-%.5f" % (lami,lamj)] = (DGij,DGij_sig)
    totline = ofile[-3].split()
    mbar = float(totline[0].strip(","))
    mbar_err = float(totline[1])
    results['DGtot'] = (mbar,mbar_err,cumtime)
        
    return results


def getNoisyWindows(energies,endtime, mythreshold=THRESHOLD):
    noisy_windows = {}
    for leg in ("bound","free"):
        noisy_windows[leg] = {}
        for stage in ("discharge","vanish"):
            noisy_windows[leg][stage] = []
            windows = list(energies[1][leg][stage][(DISCARD,endtime)].keys())
            windows.remove('DGtot')
            avg_windows = {}
            for window in windows:
                avg_windows[window] = []
                vals = []
                for run in range(1,NRUNS+1):
                    vals.append(energies[run][leg][stage][(DISCARD,endtime)][window][0])
                std = np.array(vals).std()
                avg_windows[window] = std
                if std > mythreshold:
                    noisy_windows[leg][stage].append( [window, std] )
    return noisy_windows

def calcEnergies(energies, basefolder, simprofile, startime, endtime):
    for run in range(1,NRUNS+1):
        for leg in ("bound","free"):
            for stage in ("discharge","vanish"):
                rundir = "%s/%s/run00%d/%s/output/lambda-*/simfile.dat" % (basefolder,leg,run,stage)
                simfiles = glob.glob(rundir)
                try:
                    energies[run]
                except KeyError:
                    energies[run] = {}
                try:
                    energies[run][leg]
                except KeyError:
                    energies[run][leg] = {}
                try:
                    energies[run][leg][stage] 
                except KeyError:
                    energies[run][leg][stage] = {}
                try:
                    energies[run][leg][stage][(startime,endtime)]
                except KeyError:
                    DG = doMBAR(simfiles, startime, endtime)
                    print ("MBAR on %s, start %s, end %s, DG %s " % (rundir,startime,endtime,DG))
                    energies[run][leg][stage][(startime,endtime)] = DG

def plotDGbind(energies):
    runs = list(energies.keys())
    chunks = energies[runs[0]]['bound']['discharge'].keys()
    list(chunks).sort()
    x_vals = []
    y_vals = []
    error_vals = []
    y_run1 = []
    y_run2 = []
    y_run3 = []
    y_run4 = []
    y_run5 = []
    y_runs = [y_run1, y_run2, y_run3, y_run4, y_run5]
    for chunk in chunks:
        DGbinds = []
        cumtime = 0.0  
        for run in runs:
            DG_bound_discharge = energies[run]['bound']['discharge'][chunk]['DGtot'][0]
            DG_bound_discharge_time = energies[run]['bound']['discharge'][chunk]['DGtot'][2]
            DG_bound_vanish = energies[run]['bound']['vanish'][chunk]['DGtot'][0]
            DG_bound_vanish_time = energies[run]['bound']['vanish'][chunk]['DGtot'][2]
            DG_free_discharge = energies[run]['free']['discharge'][chunk]['DGtot'][0]
            DG_free_discharge_time = energies[run]['free']['discharge'][chunk]['DGtot'][2]            
            DG_free_vanish = energies[run]['free']['vanish'][chunk]['DGtot'][0]
            DG_free_vanish_time = energies[run]['free']['vanish'][chunk]['DGtot'][2]            
            DGbind = (DG_free_discharge + DG_free_vanish) -  (DG_bound_discharge + DG_bound_vanish)
            DGbinds.append(DGbind)
            chunktime = DG_bound_discharge_time + DG_bound_vanish_time + DG_free_discharge_time + DG_free_vanish_time
            cumtime += chunktime
            y_runs[run-1].append(DGbind)
            #print (chunk,DGbind)
        DGbind_avg = np.array(DGbinds).mean()
        DGbind_ste = np.array(DGbinds).std()/math.sqrt(len(DGbinds))*1.96# 95% CI
        print (chunk[1], DGbind_avg,DGbind_ste)
        x_vals.append(cumtime)
        y_vals.append(DGbind_avg)
        error_vals.append(DGbind_ste)
    
    #ostream = open("plot.dat","w")
    #for x in range(0,len(x_vals)):
    #    line = " %8.5f " % x_vals[x]
    #    for entry in y_runs:
    #        line += " %8.5f " % entry[x]
    #    line += " %8.5f %8.5f \n" % (y_vals[x],error_vals[x])
    #    ostream.write(line)
    #ostream.close()
    ##print (x_vals)
    ## Convergence plot
    #plt.plot(x_vals,y_vals)
    #for entry in y_runs:
    #    plt.plot(x_vals,entry)
    #plt.ylabel('DG / kcal.mol-1')
    #plt.xlabel('Cumulative sampling time / ns')
    #plt.fill_between(x_vals, np.array(y_vals)-np.array(error_vals), np.array(y_vals)+np.array(error_vals), alpha=0.5, face#color='#ffa500' )
    #plt.show()

def updateSimProfile(simprofile,noisy_windows,samplingtime):
    for run in range(1,6):
        for leg in ("bound","free"):
            for stage in ("discharge","vanish"):
                simfiles = simprofile[run][leg][stage].keys()
                #print (simfiles)
                for simfile in simfiles:
                    # Ugly
                    lamval = float(os.path.split(simfile)[-2].split("/")[-1].strip("lambda-"))
                    #print (simfile)
                    for window in noisy_windows:
                        noisyleg = window[0]
                        noisystage = window[1]
                        noisylami, noisylamj = window[2].split('-')
                        noisylami = float(noisylami)
                        noisylamj = float(noisylamj)
                        if (leg == noisyleg and stage == noisystage 
                            and (abs(lamval-noisylami) < 0.01 or 
                                 abs(lamval-noisylamj) < 0.01) ):
                            print ("%s noisy, extending sampling to %s" % (simfile,samplingtime))
                            simprofile[run][leg][stage][simfile] = samplingtime
                            break
                    


#####################
basefolder = "."

energies={}
simprofile={}
# scan data to work out in which epoch we are
maxendtime = -1
for run in range(1,NRUNS+1):
    simprofile[run] = {}
    for leg in ("bound","free"):
        simprofile[run][leg] = {}
        for stage in ("discharge","vanish"):
            simprofile[run][leg][stage] = {}
            rundir = "%s/%s/run00%d/%s/output/lambda-*/simfile.dat" % (basefolder,leg,run,stage)
            simfiles = glob.glob(rundir)
            for simfile in simfiles:
                #Find how much sampling was done
                istream = open(simfile,'r')
                buffer = istream.readlines()
                istream.close()
                nsteps = int(buffer[-1].split()[0])
                endtime = (nsteps+NRGFREQ)*TIMESTEP
                if endtime > maxendtime:
                    maxendtime = endtime
                simprofile[run][leg][stage][simfile] = endtime
#print (simprofile)
print ("MAXIMUM SAMPLING TIME IS %s ns " % maxendtime)
epoch = int(maxendtime/EPOCHTIME)
print ("WE ARE IN EPOCH %s " % epoch)
#sys.exit(-1)
# Now compute energies for each sim using variable chunks
calcEnergies(energies, basefolder, simprofile, DISCARD, maxendtime)
print ("The binding free energy estimate is:")
plotDGbind(energies)

# Now prepare new submission script
templateslurm = """#!/bin/bash
#SBATCH -o adaptive-epoch-[epochnum]-job-%A.%a.out
#SBATCH -p gpu
#SBATCH -n 1
#SBATCH --gres=gpu:1
#SBATCH --time 24:00:00
#SBATCH --array=0-[endarray]

module load cuda/9.2

echo "CUDA DEVICES:" $CUDA_VISIBLE_DEVICES

lamvals=( [lamvals] )
lam=${lamvals[SLURM_ARRAY_TASK_ID]}

echo "lambda is: " $lam

cd lambda-$lam

export OPENMM_PLUGIN_DIR=/export/users/julien/sire.app/lib/plugins/

srun /export/users/julien/sire.app/bin/somd-freenrg -C ../../input/sim-extend.cfg -l $lam -p CUDA
cd ..

wait
"""

# Now work out list of noisy windows
noisy_windows = getNoisyWindows(energies,maxendtime, mythreshold=THRESHOLD)
print ("The noisy windows (std > %s kcal/mol) are:" % THRESHOLD)
for leg in noisy_windows.keys():
    for stage in noisy_windows[leg].keys():
        print ("for %s %s" % (leg,stage))
        print (noisy_windows[leg][stage])
        nwindows = len(noisy_windows[leg][stage])
        if nwindows == 0:
            continue
        lamvals = []
        for windowpair, std in noisy_windows[leg][stage]:
            winA, winB = windowpair.split("-")
            if not winA in lamvals:
                lamvals.append(winA)
            if not winB in lamvals:
                lamvals.append(winB)
        lamvals.sort()
        print (lamvals)
        lamvalstr = ""
        for lv in lamvals:
            lamvalstr += "%.3f " % (float(lv))
        for run in range(1,NRUNS+1):
            rundir = "%s/%s/run00%d/%s/" % (basefolder,leg,run,stage)
	    #sub file
            ostream = open(os.path.join(rundir,"resub-epoch-%s.sh" % (epoch+1)),'w')
            slurm = copy.deepcopy(templateslurm)
            slurm = slurm.replace("[epochnum]","%s" % (epoch+1))
            slurm = slurm.replace("[endarray]","%s" % (len(lamvals)-1))
            slurm = slurm.replace("[lamvals]",lamvalstr)
            ostream.write(slurm)
            ostream.close()
	    # config file
            ostream = open(os.path.join(rundir,"input/sim-extend.cfg"),'w')
            istream = open(os.path.join(rundir,"input/sim.cfg"),'r')
            buff = istream.readlines()
            for line in buff:
                if line.startswith("minimise"):
                    ostream.write("minimise = False\n")
                else:
                    ostream.write(line)
            ostream.close()
            istream.close()
sys.exit(-1)

