#!/usr/bin/env python3
__description__ = \
"""
Fit enrichment curves individually.
"""
__author__ = "Michael J. Harms"
__usage__ = "indivFit.py curve_pickle_file"
__date__ = "2015-01-06"

import sys, os, pickle
import numpy as np
from scipy.optimize import curve_fit
from scipy import exp

import warnings
warnings.filterwarnings('error')

def fitModel(x,decay_constant=1.0,coeff=1.0):

    return coeff*(decay_constant**x)

def doFit(x,y,decay_guess=1,coeff_guess=1.0):
    """
    """

    param, cov = curve_fit(fitModel,x,y,p0=(decay_guess,coeff_guess))

    try:
        ssres  = sum((y - fitModel(x,param[0],param[1]))**2)
        sstot  = sum((y - np.mean(y))**2)
        Rsq = 1 - ssres/sstot
    except RuntimeWarning:
        Rsq = 0.

    return param, Rsq
   
def normalizeDict(some_dict):
    """
    """

    sums = np.array([0,0,0],dtype=float)
    for key, value in some_dict.items():
        sums += np.array(value)  

    out_dict = {}
    for key in some_dict.keys():
        out_dict[key] = np.array(some_dict[key])/sums

    return out_dict

def main(argv=None):
    """
    """

    if argv == None:
        argv = sys.argv[1:]

    try:
        curve_file = argv[0]
    except IndexError:
        err = "Incorrect arguments. Usage:\n\n%s\n\n" % __usage__
        raise IndexError(err)


    f = open(curve_file,'rb')
    curve_dict = pickle.load(f)
    f.close()

    norm_curve_dict = normalizeDict(curve_dict)

    g = open("fit-output.txt",'w')
    g.write(" seq alpha beta rsq r1 r2 r3 total nr1 nr2 nr3\n")
    counter = 1
    for key, value in curve_dict.items():

        norm_counts_str = " ".join(["%e" % v for v in norm_curve_dict[key]])
        counts_str = " ".join(["%i" % v for v in value])

        y = np.array(norm_curve_dict[key])
        x = np.array(range(1,len(y)+1))

        total = sum(value)
        try:
            param, rsq = doFit(x,y,1.0,1.0/len(list(curve_dict.items())))
            g.write("%i %s %e %e %e %s %i %s\n" % (counter,key,param[0],
                                                   param[1],rsq,counts_str,
                                                   total,norm_counts_str))
        except RuntimeError:
            g.write("%i %s NA NA NA %s %i %s\n" % (counter,key,counts_str,
                                                   total,norm_counts_str))

        counter += 1

    g.close()


if __name__ == "__main__":
    main()
