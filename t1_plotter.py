import sys
import matplotlib.pyplot as pyplot
import numpy as np
import math as math
import sympy as sym
import backend as back
from scipy.optimize import curve_fit

if len(sys.argv)!=3:
	print ('Wrong number of arguments.')
	print ("Usage: " + sys.argv[0] + " <input file> " + " <output file> ")
	quit()
else:
	inputfile = sys.argv[1]    
	outfile = sys.argv[2]

output = open(outfile,"w")
infile = open(inputfile,"r")
tau = []
norm_sig = []
yer = []
found_data = False

#data_type = input('Data type to be ploted \n1.T1 \n2.T2 \n3.Rabi \n:')

while True:
	line = infile.readline()
	if '#' in line:
		break
	if found_data:
		tokens = line.split()
		tau.append(float(tokens[0]))
		norm_sig.append(float(tokens[1]))
		yer.append(float(tokens[5])/2.0)
	if '[Data]' in line:
		found_data = True

filename = input('enter a filename for the graph: ')
titlename = input('enter a title for the graph: ')

print ('Determine guess parameters that fit the formula f(x) = ae^bx + c')

a = float(input('enter the parameter a: '))
b = float(input('enter the parameter b: '))
c = float(input('enter the parameter c: '))

tau_0 = np.array(tau,dtype=float)
norm_sig_0 = np.array(norm_sig,dtype=float)
guess = (a,b,c)
popt, pcov = curve_fit(back.func,tau_0,norm_sig_0,p0=guess,maxfev=20000)
output.write("f(x) = %0.1f + %0.2fe^%gt" % (popt[2], popt[0], popt[1]))


xs = sym.Symbol('\lambda')    
tex = sym.latex(back.funcs(xs,*popt)).replace('$', '')
a2,b2,c2 = np.sqrt(np.diag(pcov))
a1,b1,c1 = popt
inv = 1.0/np.exp(1)
decay_const = 1/b1
error_const = (b2/b1)*(1/b1)
print (error_const)

percent_error = 100.0*(error_const/decay_const)
output.write("\n\ndecay constant = %.3f\n" % decay_const)
output.write("Error in decay constant = %.3f\n" % error_const)
output.write("Percent Error in decay constant = %.3f\n" % percent_error)





pyplot.figure()
pyplot.subplot(111)
pyplot.errorbar(tau,norm_sig,yerr = yer,fmt = 'x')
pyplot.plot(tau_0,back.func(tau_0,*popt), 'r-', label = ('T1 = %.3fns' % decay_const))
pyplot.title(str(titlename))
pyplot.legend(loc='best')
pyplot.xlim(xmin=0)
pyplot.xlabel('Tau time (ns)')
pyplot.ylabel('Normalised Signal')
pyplot.savefig(str(filename))

