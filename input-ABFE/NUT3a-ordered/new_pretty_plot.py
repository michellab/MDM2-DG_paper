#!/usr/bin/env python
# coding: utf-8
import sys
import numpy as np
import math
import matplotlib.pyplot as plt
       
#####################

istream = open(sys.argv[1],'r')
figfile = sys.argv[1].replace(".dat",".png")
buffer = istream.readlines()
istream.close()
nruns = len(buffer[1].split())-3
x_vals = []
y_vals = []
for y in range(0,nruns):
    y_vals.append([])
y_avg = []
error_vals = []

epoch_step = 9
epoch_counter = 1
bars = []
for line in buffer:
    if line.startswith("#"):
        continue
    elems = line.split()
    ctime = float(elems[0])
    x_vals.append(ctime)
    epoch_counter += 1
    if (epoch_counter > epoch_step):
        bars.append(ctime)
        #print ctime
        epoch_counter = 0
    for runid in range(0,nruns):
        y = float(elems[runid+1])
        y_vals[runid].append(y)
    DG_avg = float(elems[-2])
    CI = float(elems[-1])
    y_avg.append(DG_avg)
    error_vals.append(CI)
# Convergence plot
plt.plot(x_vals,y_avg)
for entry in y_vals:
    plt.plot(x_vals,entry, linestyle='dashed')
plt.ylabel('$\Delta \it{G}$ / kcal.mol$^-$$^1$')
plt.xlabel('Cumulative sampling time / ns')
plt.fill_between(x_vals, np.array(y_avg)-np.array(error_vals), np.array(y_avg)+np.array(error_vals), alpha=0.5, facecolor='#ffa500' )

for b in bars:
    plt.axvline(x=b,c='b', linestyle='dotted')
plt.savefig(figfile)
plt.show()
