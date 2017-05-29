import NWFunction_v2 as nwf
import os
# This script is part 1 of the DFT gamma optimisation
# It should be run in a directory with the .nw input file
# It will then create and run the initial jobs as laid out below

# This script creates the job environment for an nw file for gamma 0, 0.1, 0.2, 0.3, 0.4 0.5
# It also creates cationic, anionic and neutral nw input files
# It then submits them all to the cluster and works in the Scratch space

# Save environmental conditions to variables for other functions
os.chdir('..')
nwf.Environment()
gamma_values = []

master_file = nwf.MasterFileName()

# Make list of gamma values
for i in range(0,501,100):
    gamma_values.append(float(i)/1000)

# Use CompSetUp function to initialise all jobs
for gamma in gamma_values:
    nwf.CompSetUp(master_file, gamma, False)

