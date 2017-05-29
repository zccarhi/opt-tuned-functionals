import NWFunction_v2 as nwf
import numpy as np
import numpy.polynomial.polynomial as poly
import os

os.chdir('..')
# Find environmental conditions
nwf.Environment()
# Ask the user what file is being operated on
master_file = nwf.MasterFileName()
gammas = nwf.GammaEnv(master_file)
# Make list of gamma values that have not already been copied from scratch - Ones not divisible by 0.1
new_gammas = []
for item in gammas:
    if int(item * 1000) % 100 != 0:
        new_gammas.append(item)

# Copy said gammas .log files from scratch space
for gamma in new_gammas:
    filegamma = int(1000 * gamma)
    # Copy .log files from scratch space for all gamma values
    working_file = nwf.WorkingFile(master_file, filegamma, '0')
    nwf.CopyLog(working_file)
    working_file = nwf.WorkingFile(master_file, filegamma, 'm1')
    nwf.CopyLog(working_file)
    working_file = nwf.WorkingFile(master_file, filegamma, 'p1')
    nwf.CopyLog(working_file)

# Adds the gammas that are divisible by 0.1 to the new gammas list
min_gamma = min(new_gammas)
new_gammas.append(min_gamma - 0.025)
new_gammas.append(min_gamma + 0.075)
new_gammas.append(min_gamma + 0.175)
new_gammas.sort()
# Finds the HOMO vector of the gamma 0.1 file (arbitrarily chosen)
working_file = nwf.WorkingFile(master_file, 100, '0')
HOMO_vector = nwf.FindHOMO(working_file)

# Store data from JSqTable in an array
gamma_J_minimise = np.zeros((len(new_gammas) , 4))
gamma_J_minimise = nwf.JSqTable(master_file, new_gammas, HOMO_vector)
J2 = gamma_J_minimise[:,3]
index = np.argmin(J2)
gamma_min = new_gammas[index]
# Only data from the 5 points surrounding the minimum are used to fit the polyniomial to increase accuracy at the minimum
# Here we check if the point 0.05 above the min is already calculated, and if it is take from there, or move 0.025 back and plot from there
x = gamma_J_minimise[index - 2:index + 3,0]
y = gamma_J_minimise[index - 2:index + 3,3]


# Use a 4th order polynomial to fit the data
# Then find the roots of the differential to find the minimum
coefs = poly.polyfit(x, y, 4)
diff = [coefs[4] * 4, coefs[3] * 3, coefs[2] * 2, coefs[1]]
roots = np.roots(diff)
print roots

opt_gamma = 0
# Extract the roots that are real, then find the closest root to the gamma value with the lowest calculated value of j^2
for num in roots:
    if num.imag == 0:
        if abs(num.real - gamma_min) < abs(opt_gamma - gamma_min):
            opt_gamma = num.real
            
print opt_gamma
# Submit a file with the optimised gamma with gamma to 3dp.
nwf.CompSetUp(master_file, float("{0:.4f}".format(opt_gamma)), True)

# OPTIONAL - edit out if needed, also submits jobs for the gamma values immediatly surroundingthe opt_gamma to check if J^2 is minimised
opt_gamma = float("{0:.3f}".format(opt_gamma))
gamma_check = [opt_gamma - 0.002, opt_gamma - 0.001, opt_gamma, opt_gamma + 0.001, opt_gamma + 0.002]
for gamma in gamma_check:
    if int(gamma * 1000) % 25 != 0:
        nwf.CompSetUp(master_file, gamma, False)
    
        
        



