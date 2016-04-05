from libc.math cimport exp
import numpy as np
cimport numpy as np
import pandas as pd
import jellyfish as jf

import random 

cimport cython

aa_string = '*ACDEFGHIKLMNPQRSTVWYBZX'
aa_dict = dict([(a, i) for a, i in zip(aa_string, range(len(aa_string)))])

def read_matrix(file_name):
    """
    read matrix text file into dictionary for scoring.
    """

    data = open(file_name).readlines()
    seq = data[0].strip('/n/r').split()

    matrix = {}

    for line in data[1:]:
        line = line.strip('/n/r').split()
        for j in range(1, len(line)):
            b = seq[j-1]
            matrix[(line[0], b)] = exp(float(line[j]))

    return matrix


cdef class DistMatrix:
    """
    Takes a file, makes a list of sequences, and computes a distance matrix.
    Distance matrix computed using 1 - exp(probability_score) with weighted scoring.
    
    args:
            matrix: distance scoring to use in calculating distance matrix between sequences.
            
            seq (optional): if not declared, assumes file is already a list of sequences. If seq = 'no' declared then
            list of sequence is pulled out from phage data file.
            
            phage_file: the file containing the sequences to be calculated into a distance matrix.
    """
    
    cdef str _phage_file
    cdef _scoring, _matrix
    
    def __init__(self, phage_file, scoring = None):
        self._phage_file = phage_file
        self._scoring = scoring
        self._matrix = read_matrix('blosum62.txt')
        
    cpdef aminos_int(self, seq):
        
        int_seq = []
        
        for i in seq:
            for j in i:
                int_seq.append(aa_dict[j])
            
        return int_seq
        
    cpdef create_list(self):
        """
        Creates a list of sequences from a file with a k_independent greater than 1.00
        If phage file is already a list of sequences, outputs that as a list.
        """

        cdef phage_data = []

        with open(self._phage_file) as data:
            next(data)
            for line in data:
                num, seq, k_glob, theta_glob, k_ind, theta_ind = line.split()
                phage_data.append(seq) if float(k_ind) > 1.000000000e+00 else None
        return phage_data
    
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cdef double score_pairwise(self, seq1, seq2): 
        """
        calculate score between two sequences using a given scoring matrix.
        """
        cdef double score = 1.0
        cdef str i, j
        
        if self._scoring == 'hamming':
            for i, j in zip(seq1, seq2):
                score += 0.0 if i == j else 1.0
        elif self._scoring == 'damerau':
            score = jf.damerau_levenshtein_distance(seq1, seq2)
        elif self._scoring == 'weighted':
            for i, j in zip(seq1, seq2):
                score *= self._matrix[(i, j)]
            score = 1 - score   
        else:
            raise ValueError('not an option')
            
        return score
    
    @cython.wraparound(False)
    @cython.boundscheck(False)
    def calc_matrix(self):
        """
        calculates pairwise distance between sequences in a given file and returns a distance matrix.
        """
        cdef:
            np.ndarray sequences = np.array(self.create_list())
            int i, j, nrow = sequences.shape[0]
            np.ndarray[dtype=double, ndim=2] D = np.zeros((nrow, nrow))
            
            double temp

        for i in range(nrow):
            for j in range(i + 1, nrow):
                temp = self.score_pairwise(sequences[i], sequences[j])
                D[i, j] = temp
                D[j, i] = temp
                
        dist_matrix = pd.DataFrame(D, index = sequences, columns = sequences)
    
        return dist_matrix