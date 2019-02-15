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
from scipy.fftpack import fft, rfft, fftshift
import itertools


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
    output3 = open(str(titlename) + 'frequency plot data.txt','w')
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
    p = 1/100      
    for z in range(len(tau) - 1):
               if (tau[z] <= 2.0):
                    tau1.append(tau1[z] + (1/p))
                    sig.append(norm_sig[z+1])

    tau_0 = np.array(tau,dtype=float)
    for i in range(len(sig)):
        if sig[i] < 2*sig[0] :
            norm_sig1.append(sig[i])
    norm_sig_0 = np.array(sig,dtype=float)
    norm_sig_0 = np.subtract(norm_sig_0,np.average(norm_sig_0))
    #print (i)
    N = len(norm_sig_0)
    k = np.arange(N)
    T = N/101
    #print (N)
    frq = k/T
    frq = frq[range(N//2)]
    fast = fft(norm_sig_0)
    fast = fast[range(N//2)]
    freq = np.fft.fftfreq(len(norm_sig_0),p)
    print (len(frq))
    y,binEdges = np.histogram(norm_sig_0,bins=20)
    bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
    x = bincenters
    #print(str(y))
    
    
    
    gauss1 = GaussianModel(prefix='g1_')
    pars = gauss1.guess(y,x=x)
    
    pars['g1_center'].set(np.average(norm_sig_0), min=min(norm_sig_0), max=max(norm_sig_0))
    pars['g1_sigma'].set(20, min=3)
    pars['g1_amplitude'].set(300, min=10)


        


 

    
    
    mod  = gauss1
    out = mod.fit(y,pars,x=x)
    
    print (outfile)
    output.write(out.fit_report(min_correl=0.25) + '\n')
    
    xv = []
    yv = []
    for i in range(len(frq)):
        #if frq[i]<=50:
        xv.append(freq[i])
        yv.append(fast[i])
    print(len(xv))
    xv,yv = zip(*sorted(zip(xv,yv)))
    fast1 = np.array(yv)
    for i in range(len( xv)):
        output3.write(str(xv[i]) +'\t'+ str(yv[i]) +'\n')
    """
    plt.rcParams["figure.figsize"] = (5,5)
    plt.plot(y,bincenters,'.', color='b')
    plt.plot(out.best_fit,bincenters,'-',color='k',label = 'best fit')
    ax = plt.axes()
    ax.axhline(mean(sig),color = 'r',label = 'mean')
    ax.axhline((mean(sig) + math.sqrt(mean(sig))),color = 'g',label = 'shot noise upper bound')
    ax.axhline((mean(sig) - math.sqrt(mean(sig))),color = 'y',label = 'shot noise lower bound')
    plt.title(str('histogram ' + titlename))
    plt.ylabel('Signal(c/ms)')
    plt.xlabel('Frequency')
    art = []
    lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
    art.append(lgd)
    plt.savefig(outfile,additional_artists=art,
    bbox_inches="tight")
    plt.close()

    plt.rcParams["figure.figsize"] = (10,5)
    plt.plot(tau1,sig,'-',label = 'time trace',color = '#2E86C1')
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
    """
    
    plt.rcParams["figure.figsize"] = (10,5)
    plt.plot(xv,abs(fast1)**2,'-',label = 'time trace',color = '#2E86C1')
    #plt.ylim(0,5000)
    #plt.xlim(-1,50)
    plt.title(str(titlename))
    art = []
    lgd = plt.legend(loc=9, bbox_to_anchor=(0.5, -0.1), ncol=2)
    art.append(lgd)
    plt.ylabel('Signal(c/ms)^2')
    plt.xlabel('frequency(Hz)')
    plt.savefig('new_' + titlename + ' frequency.png',additional_artists=art,
    bbox_inches="tight",transparent = True)
    plt.close()
    
    del tau[:]
    del tau1[:]
    del norm_sig[:]
    del norm_sig1[:]
    del sig[:]
    del xv
    del yv

    
    i += 1
