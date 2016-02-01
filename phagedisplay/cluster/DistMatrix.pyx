from libc.math cimport exp
import numpy as np
cimport numpy as np
import pandas as pd

cimport cython

def readMatrix(fileName):
    """
    read matrix text file into dictionary for scoring.
    """

    data = open(fileName).readlines()
    seq = data[0].strip('/n/r').split()

    matrix = {}

    for line in data[1:]:
        line = line.strip('/n/r').split()
        for j in range(1, len(line)):
            b = seq[j-1]
            matrix[(line[0], b)] = int(line[j])

    return matrix

cdef class DistMatrix:
    """
    Takes a file, makes a list of sequences, and computes a distance matrix.
    Distance matrix computed using 1 - exp(probability_score).
    
    args:
            matrix: the biopython matrix to use for pairwise scoring. Can score using hamming distance 
                    if 'hamming' or weighted blosum scoring if 'weighted'.
            
            seq (optional): if not declared, assumes file is already a list of sequences. If seq = 'no' declared then
            list of sequence is pulled out from phage data file.
            
            phage_file: the file containing the sequences to be calculated into a distance matrix.
    """
    
    cdef str _phage_file, _seq
    cdef _scoring, _matrix
    
    def __init__(self, phage_file, scoring = None, seq = None):
        self._phage_file = phage_file
        self._scoring = scoring
        self._seq = seq
        _matrix = readMatrix('blosum62.txt')
        
    cpdef createSeqList(self):
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
    cdef double scorePairwise(self, seq1, seq2): 
        """
        calculate score between two sequences using a given scoring matrix.
        """
        cdef double score = 0.0
        cdef str i, j
        
        
        if self._scoring == 'hamming':
            for i, j in zip(seq1, seq2):
                score += 0.0 if i == j else 1.0
        elif self._scoring == 'weighted':
            for i, j in zip(seq1, seq2):
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

        #for i in prange(nrow, nogil = True, schedule = 'guided'):
        for i in range(nrow):
            for j in range(i + 1, nrow):
                #with gil:
                temp = self.scorePairwise(sequences[i], sequences[j])
                D[i, j] = temp
                D[j, i] = temp
                
        dist_matrix = pd.DataFrame(D, index = sequences, columns = sequences)
    
        return dist_matrix

import random 

def genSeqList(N, M, l, q):
    """
    generate dummy data.
    
    args:
        N: number of clusters
        M: differences in each sequence
        l: length of sequences
        q: minimum number of sequences in each cluster
    """
    
    aminos = ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'Y', 'W', 'E', 'D', 'K', 'R', 'H', 'S', 'T', 'C', 'M', 'N', 'Q']
    seqList = []
    
    for i in range(N):
        seq = [random.choice(aminos) for x in range(l)]
        seqList.append(''.join(seq))
        
        # randomizes number of sequences in each cluster.
        randClustNum = random.randint(q, q+20)
        for j in range(randClustNum):
            seq2 = seq[:]
            for k in range(M):
                index = random.randint(0, l-1)
                seq2[index] = random.choice(aminos)
            seqList.append(''.join(seq2))
        
    return pd.DataFrame(seqList)

def toFile(seqList, fileName):
    
    seqList.to_csv(fileName, index = False, header = False)