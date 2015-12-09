from libc.math cimport exp
import numpy as np
cimport numpy as np
from Bio.SubsMat.MatrixInfo import *
import pandas as pd
from scipy.spatial.distance import *

cimport cython
cimport openmp
from cython.parallel cimport prange

cdef class DistMatrix:
    """
    Takes a file, makes a list of sequences, and computes a distance matrix.
    Distance matrix computed using 1 - exp(probability_score).
    
    args:
            matrix: the biopython matrix to use for pairwise scoring. Can score using hamming distance 
                    if matrix = 'hamming' or a given matrix.
            
            seq (optional): if not declared, assumes file is already a list of sequences. If seq = 'no' declared then
            list of sequence is pulled out from phage data file.
            
            phage_file: the file containing the sequences to be calculated into a distance matrix.
    """
    
    cdef str _phage_file, _seq
    cdef _matrix
    
    def __init__(self, phage_file, matrix = None, seq = None):
        self._phage_file = phage_file
        self._matrix = matrix
        self._seq = seq
        
    cpdef createSeqList(self):
        """
        Creates a list of sequences from a file with a k_independent greater than 1.00
        If phage file is already a list of sequences, outputs that as a list.
        """
    
        #create 12xN array of sequences.
        
        #cdef np.ndarray[dtype=ch*, ndim=2] D = np.zeros((12, N))
        cdef phage_data = []
        
        if self._seq == 'no':
            with open(self._phage_file) as data:
                next(data)
                for line in data:
                    num, seq, k_glob, theta_glob, k_ind, theta_ind = line.split()
                    phage_data.append(seq) if float(k_ind) < 1.000000000e+00 else None
            return phage_data
        elif self._seq == 'yes':
            return [line.strip() for line in open(self._phage_file)]
        else:
            raise ValueError("invalid entry!")
    
    @cython.wraparound(False)
    @cython.boundscheck(False)
    cdef double scorePairwise(self, seq1, seq2): 
        """
        calculate score between two sequences using a given scoring matrix.
        """
        cdef double score = 0.0
        cdef str i, j
        
        
        if self._matrix == 'hamming':
            for i, j in zip(seq1, seq2):
                score += 0.0 if i == j else 1.0
        elif self._matrix:
            for i, j in zip(seq1, seq2):
                if (i, j) not in self._matrix:
                    j, i = i, j
                score += self._matrix[(i, j)]
            score = 1 - exp(score)   
        else:
            raise ValueError('not an option')

        return score
    
    @cython.wraparound(False)
    @cython.boundscheck(False)
    def calcDistMatrix(self):
        """
        calculates pairwise distance between sequences in a given file and returns a distance matrix.
        """
        cdef:
            np.ndarray sequences = np.array(self.createSeqList())
            int i, j, nrow = sequences.shape[0]
            np.ndarray[dtype=double, ndim=2] D = np.zeros((nrow, nrow))
            
            double temp

        for i in prange(nrow, nogil = True, schedule = 'guided'):
            for j in range(i + 1, nrow):
                #seq1 = np.array(list(np.ndarray[dtype=double, ndim=2] D = np.zeros((nrow, nrow))))
                #seq2 = np.array(list(sequences[j]))
                with gil:
                    temp = self.scorePairwise(sequences[i], sequences[j])
                D[i, j] = temp
                D[j, i] = temp
                
        dist_matrix = pd.DataFrame(D, index = sequences, columns = sequences)
    
        return dist_matrix