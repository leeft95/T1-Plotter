import os
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import sympy as sym
from scipy.optimize import curve_fit
from lmfit.models import GaussianModel
import matplotlib.ticker as ticker
from statistics import mean

def index_of(arrval, value):
    "return index of array *at or below* value "
    if value < min(arrval):  return 0
    return max(np.where(arrval<=value)[0])

folders = []
files = []
outname = []
i = 0
m = 0
tau = []
norm_sig = []
tau1 = []
norm_sig1 = []
sig = []
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
    outfile = 'new hist_' + os.path.splitext(outname[i])[0] + '.png'
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
        

    tau1.append(float(0.0))
    sig.append(float(norm_sig[0]))
               
    for z in range(len(tau) - 1):
               if (tau[z] <= 80.0):
                   tau1.append(tau1[z] + (1/101))
                   sig.append(norm_sig[z+1])

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
    
    
    
    gauss1 = GaussianModel(prefix='g1_')
    pars = gauss1.guess(y,x=x)
    
    pars['g1_center'].set(np.average(norm_sig_0), min=min(norm_sig_0), max=max(norm_sig_0))
    pars['g1_sigma'].set(20, min=3)
    pars['g1_amplitude'].set(300, min=10)


        


 

    
    
    mod  = gauss1
    out = mod.fit(y,pars,x=x)
    
    print (outfile)
    output.write(out.fit_report(min_correl=0.25) + '\n')
    for i in range(len(bincenters)):
        output3.write(str(bincenters[i]) + '\t' + str(y[i]) + '\t' + str(out.best_fit[i]) +'\n')
    







    plt.plot(y,bincenters,'.', color='blue')
    plt.plot(out.best_fit,bincenters,'-',color='red')
    #ax = plt.axes()
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(0.01))
    #ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0005))
    plt.title(str('histogram ' + titlename))
    plt.ylabel('Signal(c/ms)')
    plt.xlabel('Frequency')
    plt.savefig(outfile)
    plt.close()


    plt.plot(tau1,sig,'-',label = 'time trace',color = '#a2cffe')
    #plt.xlim(0,40)
    ax = plt.axes()
    ax.axhline(mean(sig),color = 'r',label = 'mean')
    ax.axhline((mean(sig) + math.sqrt(mean(sig))),color = 'g',label = 'shot noise upper bound')
    ax.axhline((mean(sig) - math.sqrt(mean(sig))),color = 'y',label = 'shot noise lower bound')
    plt.title(str(titlename))
    art = []
    lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
    art.append(lgd)
    plt.ylabel('Signal(c/ms)')
    plt.xlabel('Time(s)')
    plt.savefig('new_' + titlename + '.png',additional_artists=art,
    bbox_inches="tight")
    plt.close()

    del tau[:]
    del tau1[:]
    del norm_sig[:]
    del norm_sig1[:]
    del sig[:]
    

    
    i += 1
