#!/usr/bin/env python3

__usage__ = "yo.py logK_file dir_with_cluster_csv"

import sys, os
from math import exp

def weight_clusters(logK_file,cluster_dir):
    """
    """
    
    K_dict = {}
    f = open(logK_file,'r')
    for l in f.readlines()[1:]:
        col = l.split()

        K_dict[col[1]] = exp(float(col[2])) 
    f.close()

    clusters = [f for f in os.listdir(cluster_dir) if f[-4:] == ".csv"]

    cluster_names = []
    cluster_counts = []
    for c in clusters:

        cluster_names.append(c[:-4])
        cluster_counts.append(0)
        filename = os.path.join(cluster_dir,c)
        f = open(filename,'r')
        for l in f.readlines():
            cluster_counts[-1] += K_dict[l.strip()]
        f.close()

    total = sum(cluster_counts)
    cluster_weights = [(cluster_counts[i]/total,cluster_names[i])
                       for i in range(len(cluster_counts))]

    cluster_weights.sort(reverse=True)

    out = []
    cum_sum = 0.0
    for i in range(len(cluster_weights)):
        cum_sum += cluster_weights[i][0]
        out.append("{} {:15.5f} {:15.5f}\n".format(cluster_weights[i][1],cluster_weights[i][0],cum_sum))


    return out    
        
                      
        

def main(argv=None):
    """
    """

    if argv == None:
        argv = sys.argv[1:]

    
    try:
        logK_file = argv[0]
        cluster_dir = argv[1]
    except IndexError:
        err = "Incorrect arguments. Usage:\n\n{}\n\n".format(__usage__)
        raise IndexError(err)

    out = weight_clusters(logK_file,cluster_dir)

    return "".join(out)

if __name__ == "__main__":
    print(main())
