import sys
import matplotlib.pyplot as pyplot
import numpy as np
import math as math
import sympy as sym


def func(x,a,b,c):
	return a*np.exp(-x*b) + c

def funcs(x,a,b,c):
	return a*sym.exp(-x*b) + c
