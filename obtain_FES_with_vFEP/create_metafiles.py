import os
import sys
import matplotlib
from matplotlib import pyplot as plt
import mdtraj as md
import numpy  as np

# Readme #############################################
# 
# Salome Llabres Prat, PhD
# University of Edinburgh / 2018 
#
# This script checks the output generated each bin of the US and writes a vFEP metafile.  
# It requires the following :
# - Raw data = stored into a folder named vas 

# Variables ##############################

# Usage
usage = "This script checks the output generated each bin of the US and writes a vFEP metafile.\nIt requires the following:\n\t- Raw data = named vas.CV1.CV2.ns.prod and stored into a folder named vas\n"

# system root name
name = ''

# Description of the CVspace
# CV1
CV1min = 5
CV1max = 45+1
CV1interval=2

# CV2
CV2min = 36
CV2max = 270+1
CV2interval=8

values_x = range(CV1min, CV1max, CV1interval)
values_y = range(CV2min, CV2max, CV2interval)

# Functions ################################

def create_window_array(CV1min, CV1max, CV1interval, CV3min, CV3max, CV3interval):
    """ Creates a np.array that stores all the bins of the 2D US. It has to be set
    up manually to account for the dimensions of your system.
    (int, int, int, int, int, int ) --> (np.array)
    (CV1min, CV1max, CV1interval, CV2min, CV2max, CV2interval) --> (np.array)
    """
    
    # c stores the total number of bins. 
    c = 0

    for x in values_x:
        for y in values_y: # Add conditions if needed
            c = c + 1

    # np.array that define each bin as [target_x, target_y]
    win = np.zeros((c, 2))

    # Fill the win array. 
    count = 0

    for x in values_x:
        for y in values_y:
            win[count][0] = x # target value
            win[count][1] = y # target value
            count = count +1
    
    return win

def write_vFEP_metafile( win, maxns, flag=False ):
    """ Writes a metafile for vFEP.
    flag = Checks it there are missing va files.
    
    ( array, int, optinonal(bool) ) --> (None)
    Bin array, # ns, True/False --> None
    """
    
    # Open metafile
    RST = open(name+".vFEP.metafile."+str(maxns)+".dat","w")
    
    # Write the existing prod files located in the vas folder
    for c in range(len(win)):
        x = int(win[:,0][c])
        y = int(win[:,1][c])
        for ns in range(1,maxns+1,1):
            if os.path.isfile("vas/va."+str(x)+"."+str(y)+"."+str(ns)+".prod"):
                RST.write("vas/va."+str(x)+"."+str(y)+"."+str(ns)+".prod "+str(float(x))+" 0.5 "+str(float(y))+" 0.12185 \n")
            else:
                # Verbose output
                if flag == True:
                    print( "vas/va."+str(x)+"."+str(y)+"."+str(ns)+".prod" )
    RST.close()    
    
# Main Fuunction #############################

def main(argv):
    # get system_name
    try:
        global name 
        name = sys.argv[1] 
    except IndexError:
        print( "ERROR - No system name provided.\n\n" )
        print( usage )
        sys.exit(2)
    
    # Generate bin array
    bins = create_window_array(CV1min, CV1max, CV1interval, CV2min, CV2max, CV2interval )
    verbose=False
    
    # Write vFEP metafile
    for ns in range(1,5,1):
        if ns == 4:
           verbose =  True
        write_vFEP_metafile(bins, ns, flag=verbose)

# Main #######################################
if __name__ == "__main__":
   main(sys.argv[1:])
