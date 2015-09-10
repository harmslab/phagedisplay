__description__ = \
"""
Common functions used across the simulation sub module.
"""
__author__ = "Michael J. Harms, harmsm@gmail.com"
__date__ = "2015-04-25"

import numpy as np

def uniqueCounter(x):
    """
    Count unique elements in an array.  Like bincount but (supposedly) faster
    and gives more useful output.  
    
    Args: 
        x: one-dimensional numpy array
    Returns:
        unique: array of unique entries and (count) how often each element in
                unique was seen. 
    
    Code fragment by Eelco Hoogendoorn
    http://w3facility.org/question/numpy-frequency-counts-for-unique-values-in-an-array/
    """
        
    unique, inverse = np.unique(x, return_inverse=True)
    count = np.zeros(len(unique), np.int)
    np.add.at(count, inverse, 1)
    
    return unique, count #np.vstack(( unique, count)).T
