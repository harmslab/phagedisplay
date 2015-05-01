#!/usr/bin/env python3
__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-01-18"
__usage__ = ""

import numpy as np
import loadData, regression
import sys

import scipy, pickle

def main(argv=None):
    """
    """

    if argv == None:
        argv = sys.argv[1:]

    try:
        input_file = argv[0]
        cg_width = float(argv[1])
    except (IndexError,ValueError):
        err = "Incorrect arguments. Usage:\n\n%s\n\n" % __usage__
        raise IndexError(err)  
 
    patterns = loadData.loadData(input_file,cg_width)

    fit_pattern = np.zeros((len(patterns),patterns[0].pattern_length),dtype=float)
    degeneracy = np.zeros((len(patterns)),dtype=int)
    for i, c in enumerate(patterns):
        if len(c.degen_sequence_set) != 0:
            mean_pattern, degen = c.calcMeanPattern()
            fit_pattern[i,:] = mean_pattern
            degeneracy[i] = degen

    m = regression.fitModel(fit_pattern,degeneracy,rounds_start_at=1)
    m.runRegression()
    m.writeOut("%s.pickle" % input_file)

    theta, K, calc_values = m.returnParam()


    totals = np.zeros((3),dtype=float)
    for i, c in enumerate(patterns):
        if len(c.degen_sequence_set) != 0:
            for j in range(len(c.degen_sequence_set)):
                for k in range(len(c.degen_sequence_set[j].sequences)):
                    totals += c.degen_sequence_set[j].pattern

    for i, c in enumerate(patterns):
        if len(c.degen_sequence_set) != 0:
            for j in range(len(c.degen_sequence_set)):
                for k in range(len(c.degen_sequence_set[j].sequences)):
                    print(c.degen_sequence_set[j].sequences[k],
                          theta[i]/degeneracy[i],
                          K[i],
                          c.degen_sequence_set[j].pattern/totals,
                          calc_values[i,:]/degeneracy[i])
 

    #regression_output = pickle.load(open("%s.pickle" % input_file,"rb"))

    #print(len(patterns),len(regression_output.x))

    #for i, c in enumerate(patterns):
    #    print(regression_output.x[i],regression_output.x[i+len(patterns)])


if __name__ == "__main__":
    main()
