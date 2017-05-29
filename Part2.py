import NWFunction_v2 as nwf
import numpy as np
import os
# This script is part 2 of the DFT gamma optimisation.
# It will first copy over the .log files from the scratch space.
# It will then extract the Total DFT Energy and HOMO values from the .log files
# It will use these to find the J^2 value for each gamma
# It will then create new jobs at steps of .025 gamma from -100 to +100 around the lowest gamma value

os.chdir('..')
# Create list containing old gamma values
old_gamma = []

# Ask the user what file is being operated on
master_file = nwf.MasterFileName()
# Find environmental conditions
nwf.Environment()

# Make list of gamma values
for i in range(0,501,100):
    old_gamma.append(float(i)/1000)

for gamma in old_gamma:
    filegamma = int(1000 * gamma)
    # Copy .log files from scratch space for all gamma values
    working_file = nwf.WorkingFile(master_file, filegamma, '0')
    nwf.CopyLog(working_file)
    working_file = nwf.WorkingFile(master_file, filegamma, 'm1')
    nwf.CopyLog(working_file)
    working_file = nwf.WorkingFile(master_file, filegamma, 'p1')
    nwf.CopyLog(working_file)
                               

# Now we must find a value for vector of the HOMO

working_file = nwf.WorkingFile(master_file, 100, '0')
HOMO_vector = nwf.FindHOMO(working_file)

# Now to iterate through the .log files and find the minimum J^2
gamma_J_minimise = np.zeros((len(old_gamma), 4))
gamma_J_minimise, a, b, c = nwf.JSqTable_final(master_file, old_gamma, HOMO_vector)

print a
print b
print c
print gamma_J_minimise
# Save the J^2 values to a list, find the index of the minimum values and then find the corresponding gamma
J2 = gamma_J_minimise[:,3]
index = np.argmin(J2)
gamma_min = old_gamma[index]
new_gamma = []

# Find a new list of gamma values to perform jobs on, steps of .025 gamma avoiding repeats of those divisible by 0.1
if gamma_min - 0.1 <= 0:
    new_gamma = [0.025, 0.05, 0.075, 0.125, 0.15, 0.175]
else:
    for i in range((int(gamma_min * 1000) - 75), int(gamma_min * 1000 + 76), 25):
        if i % 100 != 0:
            new_gamma.append(float(i)/1000)

# Submit new jobs to the cluster

print new_gamma
for gamma in new_gamma:
    nwf.CompSetUp(master_file, gamma, False)
    
            

    
    
    
    
    
    
    
    

