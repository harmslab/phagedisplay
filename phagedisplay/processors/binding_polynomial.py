__description__ = \
"""
Perform a linear transformation of measured enrichment data to back out the 
enrichment factor (a binding constant, in a perfect world) for each sequence
as well as the value of the summed partition function (omega).  
"""
__author__ = "Michael J. Harms"
__date__ = "2016-06-21"

from . import BaseProcessor

import os
import numpy as np

def get_logK(prev_counts,this_counts):
    """
    """
    
    if len(prev_counts) != len(this_counts):
        err = "Count arrays must have the same number of entries.\n"
        raise ValueError(err)
        
    num_prev_zero = np.sum(prev_counts == 0)
    num_this_zero = np.sum(this_counts == 0)
    if num_prev_zero != 0 or num_this_zero != 0:
        err = "All array entries must be non-zero\n"
        raise ValueError(err)
        
    # output array [[log(K), lower_bound,upper_bound],...]
    out = np.zeros((len(prev_counts),3),dtype=float)
    
    # create index matrix and its inverse
    # [[-1,0,0],
    #  [-1,1,0],
    #  [-1,0,1]]
    # column 0 is omega (subtracted from everyone).  
    X = np.eye(len(prev_counts))
    X[:,0] = -1 
    inv_X = np.linalg.inv(X)
    
    # Determine the enrichment ratio for each set of counts
    enrich_ratio = (this_counts/np.sum(this_counts)) * (np.sum(prev_counts)/prev_counts)

    # Use linear model to solve for log(K)
    out[:,0]= np.dot(inv_X,np.log(enrich_ratio))
    
    # Now propagate uncertainty
    sigma_ratio = enrich_ratio * np.sqrt(this_counts**(-3) + prev_counts**(-3))
    K_sigma = np.sqrt(np.dot(np.abs(inv_X),sigma_ratio**2))
   
    K = np.exp(out[:,0])
    out[:,1] = np.log(K - K_sigma)
    out[:,2] = np.log(K + K_sigma)
    
    # Grab value for omega
    omega = out[0,:]
    
    # Set log(K0) to zero.  This is our reference state
    out[0,0] = 0.0
    out[0,1] = 0.0
    out[0,2] = 0.0
    
    return out, omega

def pre_treat_arrays(seq_array,prev_counts,this_counts):
    """
    Get count arrays ready.
    """
    
    # Arrays of zero counts
    prev_zero = (prev_counts == 0)
    this_zero = (this_counts == 0)
    
    # Toss data for which both previous and current counts are 0    
    at_least_one_not_zero = np.logical_not(prev_zero * this_zero)
    
    seq_array   =   seq_array[at_least_one_not_zero]
    prev_counts = prev_counts[at_least_one_not_zero]
    this_counts = this_counts[at_least_one_not_zero]
    
    # Update zero count arrays to reflect toss
    prev_zero = prev_zero[at_least_one_not_zero]
    this_zero = this_zero[at_least_one_not_zero]
    
    # Collapse everyone for which previous or current is zero into a single large bin
    either_is_zero = prev_zero + this_zero
    
    prev_zero_counts = np.sum(prev_counts[either_is_zero])
    this_zero_counts = np.sum(this_counts[either_is_zero])
    
    # Now, remove entries for wich either is zero
    neither_is_zero = np.logical_not(either_is_zero)
    
    seq_array   =   seq_array[neither_is_zero]
    prev_counts = prev_counts[neither_is_zero]
    this_counts = this_counts[neither_is_zero]
    
    # Sort array from smallest to largest number of counts
    array_order = np.argsort(prev_counts + this_counts)
    
    seq_array   =   seq_array[array_order]
    prev_counts = prev_counts[array_order]
    this_counts = this_counts[array_order]
    
    # Add final junk category to the end of the array
    if prev_zero_counts > 0 and this_zero_counts > 0:
        seq_array =   np.append(seq_array,"".join(["X" for i in range(len(seq_array[0]))]))
        prev_counts = np.append(prev_counts,prev_zero_counts)
        this_counts = np.append(this_counts,this_zero_counts)
    
    return seq_array, prev_counts, this_counts

class BindingPolynomialProcessor(BaseProcessor):

    def process(self,count_dict,
                ref_round=None,measured_round=None,
                seq_to_take=None):

        if ref_round == None:
            ref_round = 0
        if measured_round == None:
            measured_round = -1

        local_count_dict = {}
        if seq_to_take != None:
            for s in seq_to_take:
                local_count_dict[s] = count_dict[s]
        else:
            local_count_dict = count_dict

        seq_array = np.array(list(local_count_dict.keys()))
        prev_round = np.zeros(len(seq_array),dtype=int)
        this_round = np.zeros(len(seq_array),dtype=int)
        for i in range(len(seq_array)):
            prev_round[i] = local_count_dict[seq_array[i]][ref_round]
            this_round[i] = local_count_dict[seq_array[i]][measured_round]
            
        new_seq_array, new_prev, new_this = pre_treat_arrays(seq_array,
                                                             prev_round,
                                                             this_round)
        logK, omega = get_logK(new_prev,new_this)

        self._out_dict = {}
        for i in range(len(new_seq_array)):
            self._out_dict[new_seq_array[i]] = list(logK[i])

        human_readable = os.path.join(self.getProperty("expt_name"),
                                      "human-readable-summary.txt")

        f = open(human_readable,"w")
        f.write("{:>15s} {:>15s} {:>15s} {:>15s} {:>15s}\n".format(" ","seq","logK","logK_low","logK_high"))
        for i, s in enumerate(self._out_dict.keys()):
            f.write("{:15d} {:>15s} {:15.8e} {:15.8e} {:15.8e}\n".format(i,s,
                                                                        self._out_dict[s][0],
                                                                        self._out_dict[s][1],
                                                                        self._out_dict[s][2]))
        f.close()
            

    @property
    def data(self):
        """
        Return a dictionary of sequences keyed to logK.
        """

        return self._out_dict
