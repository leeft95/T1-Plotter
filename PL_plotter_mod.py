import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import math as math
import sympy as sym
from scipy.optimize import curve_fit
from lmfit.models import GaussianModel

def index_of(arrval, value):
    "return index of array *at or below* value "
    if value < min(arrval):  return 0
    return max(np.where(arrval<=value)[0])


#Here are the functions that I use (only doubleGauss is fitted - indivGauss is just
#for making the plots and takes its arguments from the fitted parameters of doubleGauss
###_________________________________________________________________###
def doubleGauss(x, a1, s1, c1, a2, s2, c2):
    "Function for the sum of two Gaussian Bells with provided arguments and parameters"
    g1 = a1*np.exp(-1*np.divide(np.square(np.divide(x-c1,s1)),2))
    g2 = a2*np.exp(-1*np.divide(np.square(np.divide(x-c2,s2)),2))
    return g1 + g2

def indivGauss(x, a1, s1, c1, a2, s2, c2):
    g1 = a1*np.exp(-1*np.divide(np.square(np.divide(x-c1,s1)),2))
    g2 = a2*np.exp(-1*np.divide(np.square(np.divide(x-c2,s2)),2))
    return g1,g2

###___________________________________________________________________###

folders = []
files = []
outname = []
i = 0
m = 0
tau = []
norm_sig = []
tau1 = []
norm_sig1 = []
yer = []
found_data = False


output = open('PL fit_stats.txt','w')
for entry in os.scandir(path = os.getcwd()):
    if entry.is_dir():
        folders.append(entry.path)
    if entry.name.startswith('time trace') and (entry.name.endswith('.png') != 1 and (entry.name.endswith('plot data.txt') !=1)):
        files.append(entry.path)
        outname.append(entry.name)

for i in range(len(files)):
    infile = open(files[i], 'r')
    outfile = 'new_' + os.path.splitext(outname[i])[0] + '.png'
    titlename = os.path.splitext(outname[i])[0]
    output.write(str(titlename) + '\n')
    output3 = open(str(titlename) + ' plot data.txt','w')
    output3.write(str(titlename) +'\n' + 'bincenters \t frequency \t bestfit \n') 
    while True:
        line = infile.readline()
        if line == '':
            found_data = False
            break
        if found_data == True:
            #print(line)
            tokens = line.split()
            tau.append(float(tokens[0]))
            norm_sig.append(float(tokens[1]))
        if 'Time' in line:
            found_data = True
        

    
    #for z in range(len(tau)):
        #tau1.append(tau[i] + (1/101))
    
    tau_0 = np.array(tau,dtype=float)
    for i in range(len(tau)):
        if norm_sig[i] < 2*norm_sig[0] :
            norm_sig1.append(norm_sig[i])
    norm_sig_0 = np.array(norm_sig1,dtype=float)    

    #print(norm_sig)

    y,binEdges = np.histogram(norm_sig_0,bins=30)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    x = bincenters
    print(str(y))
    
#Here are the initial guess values for the parameters, which are
#listed in order of their appearance in the "doubleGauss" function, but but briefly here is the heuristic:
#Guess 1/2 the max count rate for the initial amplitude of each gaussian
#Guess the variance of the data set for the variance of each gaussian
#Guess the center of the x-axis plus or minus one standard deviation for the center of each gaussian
    datavar = np.var(norm_sig_0)
    datastd = np.sqrt(datavar)
    maxval = max(y)
    p0 = [maxval/2,datavar,x[0]+((x[-1]-x[0])/2)+datastd,maxval/2,datavar,x[0]+((x[-1]-x[0])/2)-datastd]

#I DID set some bounds for the variables actually, but extremely lenient just to make sure that it doesn't blow up to
#inf and/or NaN
    fitparams, pcov = curve_fit(doubleGauss, x, y, p0, bounds = ([0,0,0,0,0,0],[maxval*100,datavar*10**2,x[-1]+1000,maxval*100,datavar*10**2,x[-1]+1000]))
    print(fitparams)
    out = doubleGauss(x, *fitparams)
    g1, g2 = indivGauss(x, *fitparams)
    NVcontrib1 = np.trapz(g1,x)
    NVcontrib2 = np.trapz(g2,x)
    NVratio = NVcontrib1/NVcontrib2

  
    print (outfile)
 #   output.write(out.fit_report(min_correl=0.25) + '\n')
    for i in range(len(bincenters)):
        output3.write(str(bincenters[i]) + '\t' + str(y[i]) + '\t' + str(out[i]) +'\n')
    output3.write('NV Ratio: '+ str(NVratio))
    







    plt.plot(y,bincenters,'.', color='blue')
    plt.plot(out,bincenters,'-',color='red')
    plt.plot(g1,bincenters,'-',color = 'green')
    plt.plot(g2,bincenters,'-',color = 'orange')
    #ax = plt.axes()
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(0.01))
    #ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0005))
    plt.title(str(titlename))
    plt.ylabel('Signal(c/ms)')
    plt.xlabel('Frequency')
    plt.savefig(outfile)
    plt.close()

    

    del tau[:]
    del tau1[:]
    del norm_sig[:]
    del norm_sig1[:] 
    

    
    i += 1
