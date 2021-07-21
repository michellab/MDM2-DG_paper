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
import pickle
### ADAPTIVE SAMPLING ROUTE
NRUNS = 5
EPOCHTIME = 5 # ns  
DISCARD = 0.5 # ns
CHUNKSIZE = 0.5 # ns
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
            time = (float(elems[0])+NRGFREQ)*TIMESTEP
            if (time < startime):
                continue
            if (time > endtime):
                break
            ostream.write(line)
        cumtime += time#(time - startime)
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


def calcEnergies(energies, basefolder, simprofile, startime, endtime):
    for run in range(1,NRUNS+1):
        for leg in ("bound","free"):
            for stage in ("discharge","vanish"):
                rundir = "%s/%s/run00%d/%s/output/lambda-*/simfile.dat" % (basefolder,leg,run,stage)
                simfiles = glob.glob(rundir)
                # use analyse_freenrg to get free energy of chunks
                nchunks = int(endtime/CHUNKSIZE)
                for chunk in range(1,nchunks):
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
                    startime = DISCARD
                    chunkendtime = (chunk+1)*CHUNKSIZE
                    chunkendtime = round(chunkendtime,2)
                    try:
                        DG = energies[run][leg][stage][(startime,chunkendtime)]
                        print ("USING cached MBAR results for %s, start %s, end %s, DG %s " % (rundir,startime,chunkendtime,DG))
                    except KeyError:
                        # Check if all windows have a maxendtime < endtime, if so use previous endtime result
                        maxsimfiletime = -1
                        for simfile in simfiles:
                            simfiletime = simprofile[run][leg][stage][simfile]
                            simfiletime = round(simfiletime,2)
                            if maxsimfiletime < simfiletime:
                                maxsimfiletime = simfiletime
                        if maxsimfiletime < chunkendtime:
                            DG = energies[run][leg][stage][(startime,maxsimfiletime)]
                            print ("No additional data to process over this time interval (%s-%s)" % (startime, chunkendtime))
                            energies[run][leg][stage][(startime,chunkendtime)] = DG
                            print ("reusing cached MBAR results %s, start %s, end %s, DG %s " % (rundir,startime,chunkendtime,DG))
                        else:
                            # Otherwise do MBAR from scratch
                            DG = doMBAR(simfiles, startime, chunkendtime)
                            print ("DOING NEW MBAR on %s, start %s, end %s, DG %s " % (rundir,startime,chunkendtime,DG))
                            energies[run][leg][stage][(startime,chunkendtime)] = DG
                        

def plotDGbind(energies, nepochs=0):
    runs = list(energies.keys())
    runs.sort()
    chunks = list(energies[runs[0]]['bound']['discharge'].keys())
    chunks.sort()
    x_vals = []
    y_vals = []
    error_vals = []
    y_runs = []
    for x in range(0,NRUNS):
        y_runs.append([])
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
    
    ostream = open("convergence-epoch-%s.dat" % nepochs,"w")
    header = "# cumulative_time DGrunX... avgDG 95CI\n"
    ostream.write(header)
    for x in range(0,len(x_vals)):
        line = " %8.5f " % x_vals[x]
        for entry in y_runs:
            line += " %8.5f " % entry[x]
        line += " %8.5f %8.5f \n" % (y_vals[x],error_vals[x])
        ostream.write(line)
    ostream.close()

    ##print (x_vals)
    ## Convergence plot
    #plt.plot(x_vals,y_vals)
    #for entry in y_runs:
    #    plt.plot(x_vals,entry)
    #plt.ylabel('DG / kcal.mol-1')
    #plt.xlabel('Cumulative sampling time / ns')
    #plt.fill_between(x_vals, np.array(y_vals)-np.array(error_vals), np.array(y_vals)+np.array(error_vals), alpha=0.5, face#color='#ffa500' )
    #plt.show()

         
#####################
basefolder = "./"
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
# Load cached energies (if any)
if os.path.exists("energies.pickle"):
    energies = pickle.load(open("energies.pickle","rb"))
else:
    energies = {}
#sys.exit(-1)
# Now compute energies for each sim using variable chunks
calcEnergies(energies, basefolder, simprofile, DISCARD, maxendtime)
plotDGbind(energies, nepochs=epoch)
# Save energies processed 
dumpme = open("energies.pickle","wb")
pickle.dump(energies, dumpme)
dumpme.close()

