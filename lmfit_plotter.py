import os
import sys
from pprint import pprint
#import matplotlib.pyplot as pyplot
import numpy as np
import math as math
import sympy as sym
import backend as back
from scipy.optimize import curve_fit
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
#output.write("Particle\t\t\t\tT1\t\t\tError\t\t\tPower\n")
i = 0
m = 0
tau = []
norm_sig = []
tau1 = []
norm_sig1 = []
yer = []
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
    params.add('b', value= 1e-5, min=1e-6, max=1e-3)
    params.add('c', value= 0.7, min=0, max=1) 
    params.add('d', value= 0.9, min=-1, max=1)
    print(titlename)
    result = minimize(func,params,args=(tau_0,norm_sig_0),scale_covar=True)
    #print(str(result.ier()))
    output.write(titlename + '\n')
    result.params.pretty_print(oneline=False, colwidth=8, precision=4, fmt='g', columns=['value','stderr'])
    final = norm_sig_0 + result.residual
    #print (fit_report(result,show_correl=False))
    output.write(fit_report(result,show_correl=False))
    #output.write(y)
    output.write('\n\n')
    
    
    try:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.subplot(111)
        plt.plot(tau_0, norm_sig_0, 'k+')
        plt.errorbar(tau,norm_sig,yerr = yer,fmt = 'x')
        plt.plot(tau_0, final, 'r')
        plt.title(str(titlename))
        plt.legend(loc='best')
        plt.xlim(xmin=0)
        plt.xlabel('Tau time (ns)')
        plt.ylabel('Normalised Signal')
        plt.savefig(str(filename))
        plt.close()
        
    except:
        pass
    
    """
    
    
    print('\n')

    
    
    
    guess = (a,b,c,d)
    popt, pcov = curve_fit(back.func,tau_0,norm_sig_0,p0=guess,maxfev=200000)
    #output.write("f(x) = %0.1f + %0.2fe^%gt" % (popt[2], popt[0], popt[1]))


    #xs = sym.Symbol('\lambda')    
    #tex = sym.latex(back.funcs(xs,*popt)).replace('$', '')
    a2,b2,c2,d2 = np.sqrt(np.diag(pcov))
    a1,b1,c1,d1 = popt
    #inv = 1.0/np.exp(1)
    print (b1)
    decay_const = 1/b1
    power = c1
    error_const = (b2/b1)*(1/b1)
    #print (error_const)

    percent_error = 100.0*(error_const/decay_const)
    output.write(titlename + "\t\t\t%.3f\t\t\t%.3f\t\t\t%.3f\n" % (decay_const,error_const,power))
    #output.write("Error in decay constant = %.3f\n" % error_const)
    #output.write("Percent Error in decay constant = %.3f\n" % percent_error)

    pyplot.figure()
    pyplot.subplot(111)
    pyplot.errorbar(tau,norm_sig,yerr = yer,fmt = 'x')
    pyplot.plot(tau_0,back.func(tau_0,*popt), 'r-', label = ('T1 = %.3fns' % decay_const))
    #pyplot.plot(tau,y_opt, 'r-', label = ('T1 = %.3fns' % decay_const))
    pyplot.title(str(titlename))
    pyplot.legend(loc='best')
    pyplot.xlim(xmin=0)
    pyplot.xlabel('Tau time (ns)')
    pyplot.ylabel('Normalised Signal')
    """
    #Fix code to histogram data
    """
    #hist, bins=np.histogram(rdf, bins=100,density=True)
    #bincenters = 0.5*(bins[1:]+bins[:-1])
    pyplot.savefig(str(filename))
    pyplot.close()
    """
    del tau[:]
    del norm_sig[:]
    del norm_sig1[:]
    del yer[:]
     
    


    i += 1
    
