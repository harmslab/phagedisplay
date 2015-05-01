#!/usr/bin/env python3
__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-01-17"
__usage__ = "clusterer.py read_pickle_file cg_width"

import sys, pickle, scipy, copy
import numpy as np
from math import sqrt

class DegenerateSequences:
    """
    Class that holds onto a set of sequences that share the exact same count
    pattern.
    """

    def __init__(self,pattern,sequence):
        """
        Initialize with a pattern and first sequence. 
        """
        self.pattern = tuple(pattern)
        self.sequences = [sequence]

    def append(self,sequence):
        """
        Add a new sequence to the list.
        """
        self.sequences.append(sequence)
        
   
class SequenceCluster:
    """
    Class that holds a set of DegenerateSequences instances that were placed 
    into a cluster because their patterns were similar.  
    """

    def __init__(self):
        """
        Create a list to hold all of the DegenerateSequences instances.
        """

        self.degen_sequence_set = []
        self.patterns_loaded = False
        self.pattern_length = -1

    def addSequenceSet(self,s):
        """
        Add a DegenerateSequences instance to the list.
        """
    

        if not self.patterns_loaded:
            self.patterns_loaded = True
            self.pattern_length = len(s.pattern)
        else:
            if len(s.pattern) != self.pattern_length:
                err = "Patterns must all have same length."
                raise ValueError(err)

        self.degen_sequence_set.append(copy.copy(s))

    def calcMeanPattern(self):
        """
        """

        pattern_list = []
        for s in self.degen_sequence_set:
            pattern_list.append(s.pattern)

        unique = list(dict([(p,()) for p in pattern_list]).keys())

        self.mean_pattern = np.zeros((len(self.degen_sequence_set[0].pattern)),
                                     dtype=float)
        self.degeneracy = 0
        for d in self.degen_sequence_set:
            pattern = np.array(d.pattern)
            repeats = len(d.sequences)
            
            self.mean_pattern = self.mean_pattern + pattern*repeats
            self.degeneracy += repeats
            
        self.mean_pattern = self.mean_pattern/self.degeneracy

        return self.mean_pattern, self.degeneracy    

    def print(self):
        print("SUMM:",self.mean_pattern,self.degeneracy,len(self.degen_sequence_set))
        for x in self.degen_sequence_set:
            for y in x.sequences:
                print("----",y,x.pattern)

class CoarseGrainer:
    """
    Coarse-grain a set of raw counts into bins.  Bin size is based on identifying
    non-overlapping bins based on poisson counting error.  In other words:

         (error)  (error)
    bin1--------||--------bin2---------|

    Bin size increases as counts increase because, for a poisson counting process,
    the absolute error increases with inreased magnitude.  Error is estimated by
    sqrt(N). 
    """   

    def __init__(self,maximum_value,cg_width=1):
        """
        Initialize self.breaks to hold a set of breaks with which to classify 
        inputs.  Maximum value sets the highest value that can be coarse-
        grained.  Anything above this will be set to the maximum bin.  Scalar
        is multipled by sqrt(value) to set width of coarse graining.  
        """

        v = np.array(range(2,maximum_value+1),dtype=int)
        v_err_floor = v - cg_width*np.sqrt(v)

        self.breaks = [0,1]
        current_ceiling = self.breaks[-1] + cg_width*sqrt(self.breaks[-1])
        for i in range(len(v)):
            if current_ceiling < v_err_floor[i]:
                self.breaks.append(v[i])
                current_ceiling = v[i] + cg_width*sqrt(v[i])

    def coarseGrain(self,v):
        """
        Walk through the breaks in self.breaks and figure out where each index
        in a list lands in the breaks list.  Pretty hacked...
        """

        # Bins holds bin assinment for each position in v; v_indexes holds the
        # indexes in v that have *not* been already assigned.  0 is always 0 in
        # our scheme, so these start as assigned
        bins = [0 for j in range(len(v))]
        v_indexes = [x for x in range(len(v)) if v[x] != 0]
       
        # Go through the breaks... 
        for i in range(1,len(self.breaks)):

            # Go through the sites in v that are not yet binned
            for j in v_indexes:
                
                # If we're in the bin, call it
                if v[j] > self.breaks[i-1] and v[j] <= self.breaks[i]:
                    bins[j] = i
                    v_indexes.remove(j)

            if len(v_indexes) == 0:
                break

        # If we didn't call everything, they are higher than max; record them as
        # being in the highest bin
        if len(v_indexes) != 0:
            for j in v_indexes:
                bins[j] = len(self.breaks)-1

        return bins           

    def bins2counts(self,bins):
        """
        Convert a set of bins spit out by coarseGrain back into a set 
        of putative counts (coarse-grained, obviously). 
        """

        return [self.breaks[b] for b in bins] 


def loadData(input_file,cg_width=1.0,minimum_times_seen=4):
    """
    """
    
    # Load file
    print("Loading %s" % input_file)
    count_dict = pickle.load(open(input_file,'rb'))
    all_sequences = list(count_dict.keys())
    num_rounds = len(count_dict[all_sequences[0]])

    # Create a list, (degen_sequence_list) that has a list of all sequences that
    # share exactly the same count pattern.  Each unique pattern/sequence set is
    # stored in an instance of the DegenerateSequences class.
    degen_sequence_dict = {}
    for s in all_sequences:
        if sum(count_dict[s]) < minimum_times_seen:
            continue

        try:
            degen_sequence_dict[tuple(count_dict[s])].append(s)
        except KeyError:
            pattern = tuple(count_dict[s])
            degen_sequence_dict[pattern] = DegenerateSequences(pattern,s)

    degen_sequence_list = []
    for k in degen_sequence_dict.keys():
        degen_sequence_list.append(copy.copy(degen_sequence_dict[k]))

    num_unique_patterns = len(degen_sequence_list)
    print("Found %i unique patterns for %i clones" % (num_unique_patterns,
                                                      len(all_sequences)))

    # Create a numpy array of all unique patterns
    unique_patterns = np.zeros((num_unique_patterns,num_rounds),dtype=int)
    for i in range(num_unique_patterns):
        unique_patterns[i,:] = np.array(degen_sequence_list[i].pattern)

    # Coarse grain these unique patterns into similar patterns
    cg = CoarseGrainer(np.max(unique_patterns),cg_width)
    similar_pattern_dict = {}
    for i in range(num_unique_patterns):

        key = tuple(cg.coarseGrain(unique_patterns[i,:]))
        try:
            similar_pattern_dict[key].append(degen_sequence_list[i])
        except KeyError:
            similar_pattern_dict[key] = [degen_sequence_list[i]]

    # Create a final set of clusters (made of SequenceCluster instance)
    final_clusters = [SequenceCluster() for i in range(len(similar_pattern_dict))]
    for i, key in enumerate(similar_pattern_dict.keys()):
        for cluster in similar_pattern_dict[key]:
            final_clusters[i].addSequenceSet(cluster)

    print("Collapsed to %i similar patterns." % (len(final_clusters)))

    return final_clusters
    
