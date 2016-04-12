#!/usr/bin/env python3
__description__ = \
"""
Converts a blosum62-style matrix into a pretty-print python nested dictionary.
"""
__author__ = "Michael J. Harms"
__date__ = "2016-04-11"
__usage__ = "convert-blosum-to-python.py blosum_file matrix_name"

import sys

def read_weight_file(weight_file):
    """
    Read a weight/distance matrix file into a numpy array.  
    """

    data = open(weight_file).readlines()

    # Grab alphabet bits from top line
    seq = data[0].strip('/n/r').split()
   
    # Go through each line, populating weight matrix. 
    weight_matrix = dict([(s,dict([(s,0) for s in seq])) for s in seq])

    for line in data[1:]:
        line = line.strip('/n/r').split()

        a = line[0]
        for j in range(1, len(line)):
            b = seq[j-1]
            weight_matrix[a][b] = float(line[j])

    return weight_matrix

def main(argv=None):
    """
    Main function.
    """

    if argv == None:
        argv = sys.argv[1:]

    try:
        input_file = argv[0]
        matrix_name = argv[1]
    except IndexError:
        err = "Incorrect arguments. Usage:\n\n{}\n\n".format(__usage__)
        raise IndexError(err)

    x = read_weight_file(input_file)

    # Print out in pretty-looking form
    keys = list(x.keys())
    keys.sort()
    print(matrix_name,"= { \\")
    for k1 in keys:
 
        print("    '{:}'".format(k1),": {",end="")
        for i, k2 in enumerate(keys):
            if i != 0 and i % 5 == 0:
                print("")
                print("           ",end="")       
 
            print("'{:s}' : {:-4.1f}, ".format(k2,x[k1][k2]),end="")
        print("},\n")
 
    print("}")

if __name__ == "__main__":
    main()
