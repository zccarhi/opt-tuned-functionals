import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import NWFunction_v2 as nwf
import numpy.polynomial.polynomial as poly

os.chdir('..')

nwf.Environment()
master_file = nwf.MasterFileName()

allgammas = nwf.GammaEnv(master_file)
print allgammas

filegammas = ['OPT']
for gamma in allgammas:
    filegammas.append(int(1000*gamma))

for filegamma in filegammas:
    working_file = nwf.WorkingFile(master_file, filegamma, '0')
    nwf.CopyLog(working_file)
    working_file = nwf.WorkingFile(master_file, filegamma, 'm1')
    nwf.CopyLog(working_file)
    working_file = nwf.WorkingFile(master_file, filegamma, 'p1')
    nwf.CopyLog(working_file)

working_file = nwf.WorkingFile(master_file, 100, '0')
HOMO_vector = nwf.FindHOMO(working_file)

filegammas.remove('OPT')
polyfitgammas = []

for a in filegammas:
    if a % 100 != 0 and a % 25 == 0:
        polyfitgammas.append(float("{0:.3f}".format(float(a)/1000)))
min_gamma = min(polyfitgammas)
polyfitgammas.append(float("{0:.3f}".format(min_gamma - 0.025)))
polyfitgammas.append(float("{0:.3f}".format(min_gamma + 0.075)))
polyfitgammas.append(float("{0:.3f}".format(min_gamma + 0.175)))
polyfitgammas.sort()
print polyfitgammas

polyJ, polyminus, polyneutral, polyplus = nwf.JSqTable_final(master_file, polyfitgammas, HOMO_vector)

print polyJ

J2 = polyJ[:,3] 
index = np.argmin(J2)
x = polyJ[index - 2: index + 3, 0]
y = polyJ[index - 2: index + 3, 3]
gamma_min = polyJ[index,0]

coefs = poly.polyfit(x, y, 4)
diff = [coefs[4] * 4, coefs[3] * 3, coefs[2] * 2, coefs[1]]
roots = np.roots(diff)
print coefs
print roots
opt_gamma = 0
# Extract the roots that are real, then find the closest root to the gamma value with the lowest calculated value of j^2
for num in roots:
    if num.imag == 0:
        if abs(num.real - gamma_min) < abs(opt_gamma - gamma_min):
            opt_gamma = num.real

gamma_check = [float("{0:.3f}".format(opt_gamma)) - 0.002, float("{0:.3f}".format(opt_gamma)) - 0.001, float("{0:.3f}".format(opt_gamma)), float("{0:.3f}".format(opt_gamma)) + 0.001, float("{0:.3f}".format(opt_gamma)) + 0.002]

opt_gamma = float("{0:.4f}".format(opt_gamma))
print opt_gamma
print gamma_check
gamma_check.append(opt_gamma)

gamma_J_minimise, minus, neutral, plus = nwf.JSqTable_final(master_file, gamma_check, HOMO_vector)
print minus
print neutral
print plus
print gamma_J_minimise
HOMO_LUMO = neutral[-1,3]-neutral[-1,2]
os.mkdir('%s_Analysis' %master_file)

f = open('%s_Analysis/%s_analysis' %(master_file, master_file), 'w')

f.write('Optimised Gamma = %s\n' %gamma_J_minimise[-1,0])
f.write('HOMO Eigenvalue = %s\n' %neutral[-1,2])
f.write('Ionisation Energy = %s\n' %gamma_J_minimise[-1,1])
f.write('LUMO Eigenvalue = %s\n' %neutral[-1,3])
f.write('HOMO-LUMO Gap = %s\n\n\n' %HOMO_LUMO)
f.write('J^2 Values of points nearby\n')
f.write('Gamma %s\n' %gamma_J_minimise[:,0])
f.write('J^2   %s' %gamma_J_minimise[:,3])

f.close()

x_new = np.arange(gamma_check[0] - 0.05, gamma_check[-1] + 0.05, 0.002)
ffit = poly.polyval(x_new, coefs)
plt.plot(x_new, ffit, 'r-')
plt.plot(x, y, 'ro')
plt.plot(gamma_J_minimise[-1,0], gamma_J_minimise[-1,3], 'bo')

plt.savefig("%s_Analysis/%s.pdf" %(master_file, master_file))

"""
gamma_J_minimise = np.zeros((len(gammas) , 4))
gamma_J_minimise, gamma_table_minus, gamma_table_neutral, gamma_table_plus = nwf.JSqTable_final(master_file, gammas, HOMO_vector)

print gamma_J_minimise


25Jmin,a,b,c = nwf.JSqTable_final(master_file, fitting_gammas, HOMO_vector)


J2 = 25Jmin[:,3]
index = np.argmin(J2)
gamma_min = new_gammas[index]
    
x = 25Jmin[index - 2:index + 3,0]
y = 25Jmin[index - 2:index + 3,3]

x_opt = gamma_J_minimise[-1,0]
y_opt = gamma_J_minimise[-1,3]

coefs = poly.polyfit(x, y, 4)
diff = [coefs[4] * 4, coefs[3] * 3, coefs[2] * 2, coefs[1]]
roots = np.roots(diff)
opt_gamma = 0
# Extract the roots that are real, then find the closest root to the gamma value with the lowest calculated value of j^2
for num in roots:
    if num.imag == 0:
        if abs(num.real - gamma_min) < abs(opt_gamma - gamma_min):
            opt_gamma = num.real
            
print opt_gamma

x_new = np.arange(25Jmin[index-2] - 0.5, 25Jmin[index+1] + 0.5, 0.002)
ffit = poly.polyval(x_new, coefs)
plt.plot(x_new, ffit, 'r-')
plt.plot(x, y, 'ro')
plt.plot(x_opt, y_opt, 'bo')plt.savefig('%s.pdf' %master_file)

plt.savefig("%s.pdf" %master_file)
"""
