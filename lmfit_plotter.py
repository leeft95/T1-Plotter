import os
import sys
from pprint import pprint
import numpy as np
import math as math
import sympy as sym
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from lmfit import *

def func(params,x,data):
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value
    d = params['d'].value
    model = a*np.exp(-(b*x)**c) + d
    return model - data
folders = []
files = []
outname = []
output = open('T1_data.txt','w')
output2 = open('T1_Fit_Stats.txt','w')
output.write("Particle\tT1(\u00B5s)\tError\tPerc Error\tPower\tError\tAmplitude\tError\n")
i = 0
m = 0
tau = []
norm_sig = []
tau1 = []
norm_sig1 = []
yer = []
dec = []
dec_1 = []
amp = []
amp_1 = []
found_data = False
#print ('Determine guess parameters that fit the formula f(x) = ae^(bx)^d + c')
#b = float(input('enter the parameter b: '))
for entry in os.scandir(path = os.getcwd()):
    if entry.is_dir():
        folders.append(entry.path)
    if entry.name.startswith('NVspin') and (entry.name.endswith('raw.txt') != 1) and (entry.name.endswith('.png') !=1) and (entry.name.startswith('NVspin_Rabi') !=1):
        files.append(entry.path)
        outname.append(entry.name)

for i in range(len(files)):
    infile = open(files[i], 'r')
    outfile = 'new_' + os.path.splitext(outname[i])[0] + '.png'

    while True:
        line = infile.readline()
        if found_data:
            #print(line)
            tokens = line.split()
            tau.append(float(tokens[0]))
            norm_sig.append(float(tokens[1]))
            yer.append(float(tokens[5])/2.0)
            m +=1
        if '[Data]' in line:
            found_data = True
        if m == 11:
            found_data = False
            m = 0
            break
       
    filename = outfile
    titlename =  os.path.splitext(outname[i])[0]
    tau_0 = np.array(tau,dtype=float)
    norm_sig_0 = np.array(norm_sig,dtype=float)

    #a,c,d = 0.1,0.7,0.8

    params = Parameters()
    params.add('a', value= 0.1, min=-1, max=1) 
    params.add('b', value= 1e-5, min=1e-8, max=1e-3)
    params.add('c', value= 0.5, min=0.4, max=0.8) 
    params.add('d', value= 0.9, min=0, max=1)
    print(titlename)
    result = minimize(func,params,args=(tau_0,norm_sig_0),scale_covar=True)
    #print(result.params['b'].value)
    

    decay_const = (1/result.params['b'].value)/1000
    dec.append(decay_const)
    amp.append(result.params['d'].value)
    error_const = (((result.params['b'].stderr/result.params['b'].value)*(1/result.params['b'].value)))/1000
    print(error_const)
    perc_error = (error_const/decay_const)*100
    output.write(titlename + '\t' + str(decay_const)+ '\t'+ str(error_const) + '\t' + str(perc_error) + '\t' +str(result.params['c'].value) +'\t'+ str(result.params['c'].stderr) + '\t' + str(result.params['d'].value) + '\t' + str(result.params['d'].stderr) + '\n')
    #print(str(result.ier()))
    #result.params.pretty_print(oneline=False, colwidth=8, precision=4, fmt='g', columns=['value','stderr'])
    final = norm_sig_0 + result.residual
    #print (fit_report(result,show_correl=False))
    output2.write(titlename + '\n')
    output2.write(fit_report(result,show_correl=False)+'\n\n')
    #output.write(y)
    
    
    
    plt.figure()
    plt.subplot(111)
    plt.plot(tau_0, norm_sig_0, 'k+')
    plt.errorbar(tau,norm_sig,yerr = yer,fmt = 'x')
    plt.plot(tau_0, final, 'r',label=('T1 = %.3f\u00B5s' % decay_const))
    plt.title(str(titlename))
    plt.legend(loc='best')
    plt.xlim(xmin=0)
    plt.xlabel('Tau time (ns)')
    plt.ylabel('Normalised Signal')
    plt.savefig(str(filename))
    plt.close()
    
    del tau[:]
    del norm_sig[:]
    del norm_sig1[:]
    del yer[:]

     
    


    i += 1

for i in range(len(dec)):
    if dec[i] > 0.002 and dec[i] < 0.7:
        dec_1.append(dec[i])

for i in range(len(amp)):
    if amp[i] > 0.6 and amp[i] < 0.9:
        amp_1.append(amp[i])
        
dec_2 = np.array(dec_1,dtype=float)
print(dec_2)
amp_2 = np.array(amp_1,dtype=float)
plt.figure()
plt.hist(dec_2,edgecolor='black', linewidth=1,bins = 10)

y,binEdges = np.histogram(dec_2,bins=100)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
#plt.plot(bincenters,y,'-', color='black')
ax1 = plt.axes()
#ax1.xaxis.set_major_locator(ticker.MultipleLocator(100))
#ax1.xaxis.set_minor_locator(ticker.MultipleLocator(20))
#plt.xlim(xmin=0)
plt.title('T1 Spread')
plt.xlabel('T1(\u00B5s)')
plt.ylabel('Frequency')
plt.savefig('T1 Histogram')
plt.close()
plt.figure()
x,binEdges_0 = np.histogram(amp_2,bins=20)
bincenters_0 = 0.5*(binEdges_0[1:]+binEdges_0[:-1])
#plt.plot(bincenters_0,x,'--',marker='x', color='black')
plt.hist(amp_2,edgecolor='blue', linewidth=1,bins = 5)
#ax = plt.axes()
#ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.0005))
plt.title('Amplitude Spread')
plt.xlabel('Amplitude')
plt.ylabel('Frequency')
plt.savefig('Amplitude Histogram')
plt.close()
