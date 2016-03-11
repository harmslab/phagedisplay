from libc.math cimport exp
import numpy as np
cimport numpy as np
import pandas as pd
import jellyfish as jf

import random 

cimport cython

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
    
    cdef str _phage_file, _seq
    cdef _scoring, _matrix, _amino_num
    
    def __init__(self, phage_file, scoring = None, seq = None):
        self._phage_file = phage_file
        self._scoring = scoring
        self._seq = seq
        self._matrix = read_matrix('blosum62.txt')
        self._amino_num = {'G' : 0, 'A' : 1, 'V': 2, 'L' : 3, 'I' : 4, 'P' : 5, 'F' : 6, 'Y' : 7, 'W' : 8, 
                           'E' : 9, 'D' : 10, 'K' : 11, 'R' : 12, 'H' : 13, 'S' : 14, 'T' : 15, 'C' : 16, 
                           'M' : 17, 'N' : 18, 'Q' : 19}
        
    
    cpdef aminos_int(self, seq):
        
        int_seq = []
        
        for i in seq:
            for j in i:
                int_seq.append(self._amino_num[j])
            
        return int_seq
        
    cpdef create_list(self):
        """
        Creates a list of sequences from a file with a k_independent greater than 1.00
        If phage file is already a list of sequences, outputs that as a list.
        """

        cdef phage_data = []
        
        if self._seq == 'no':
            with open(self._phage_file) as data:
                next(data)
                for line in data:
                    num, seq, k_glob, theta_glob, k_ind, theta_ind = line.split()
                    phage_data.append(seq) if float(k_ind) > 1.000000000e+00 else None
            return phage_data
        elif self._seq == 'yes':
            return [line.strip() for line in open(self._phage_file)]
        else:
            raise ValueError("invalid entry!")
    
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
        if self._scoring == 'damerau':
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

def gen_list(l, file_name):
    """
    generate dummy data.
    
    args:
        N: number of clusters
        M: differences in each sequence
        l: length of sequences
        q: max number of sequences in each cluster
    """
    
    aminos = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'E', 'D', 'K', 'R', 'H', 'S', 'T', 'C', 'M', 'N', 'Q']
    seq_list = []
    
    N = random.randint(50, 150)
    M = random.randint(2, 7)
    
    file_save = {'# clusters': N, '# of differences' : M}
    pickle.dump(file_save, open(file_name, 'wb'))
    
    for i in range(N):
        seq = [random.choice(aminos) for x in range(l)]
        seq_list.append(''.join(seq))
        
        # randomizes number of sequences in each cluster.
        rand_clust_num = random.randint(5, 250)
        for j in range(rand_clust_num):
            seq2 = seq[:]
            for k in range(M):
                index = random.randint(0, l-1)
                seq2[index] = random.choice(aminos)
            seq_list.append(''.join(seq2))
        
    return pd.DataFrame(seq_list)

def to_file(seq_list, file_name):
    
    seq_list.to_csv(file_name, index = False, header = False)