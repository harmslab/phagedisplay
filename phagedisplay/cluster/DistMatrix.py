from Bio.SubsMat.MatrixInfo import *
from sklearn.preprocessing import *
from skbio import DistanceMatrix
import numpy as np
import pandas as pd
import os

class DistMatrix():
    """
    Takes a file, makes a list of sequences, and computes a distance matrix.
    Distance matrix computed using 1 - exp(probability_score).
    
    args
            matrix: the biopython matrix to use for pairwise scoring. Can score using hamming distance 
                    if matrix = 'hamming' or a given matrix.
            
            seq (optional): if not declared, assumes file is already a list of sequences. If seq = 'no' declared then
            list of sequence is pulled out from phage data file.
            
            phage_file: the file containing the sequences to be calculated into a distance matrix.
    """
    
    def __init__(self, phage_file, matrix = None, seq = None):
        self._phage_file = phage_file
        self._matrix = matrix
        self._seq = seq
        
    def createSeqList(self):
        """
        Creates a list of sequences from a file with a k_independent greater than 1.00
        If phage file is already a list of sequences, outputs that as a list.
        """
    
        phage_data = []
        
        if self._seq == 'no':
            with open(self._phage_file) as data:
                next(data)
                for line in data:
                    num, seq, k_glob, theta_glob, k_ind, theta_ind = line.split()
                    phage_data.append(seq) if float(k_ind) > 1.000000000e+00 else None
            return phage_data
        else:
            return [line.strip() for line in open(self._phage_file)]
    
    def pairCheck(self, i, j):
        """
        check if pair (i, j) is in biopython matrix, if not returns reversed pair.
        """
        
        if (i, j) not in self._matrix:
            return self._matrix[(tuple(reversed((i,j))))]
        else:   
            return self._matrix[(i, j)]
    
    def scorePairwise(self, seq1, seq2): 
        """
        calculate score between two sequences using a given scoring matrix.
        """
        
        score = 0
        
        if self._matrix == 'hamming':
            for i, j in zip(seq1, seq2):
                score += 0 if i == j else 1
        else:
            for i, j in zip(seq1, seq2):
                score += self.pairCheck(i, j)
            score = 1 - np.exp(score)
                

        return score

    
    def calcDistMatrix(self):
        """
        calculates pairwise distance between sequences in a given file and returns a distance matrix.
        """
    
        sequences = np.array(self.createSeqList())
        dist_func = lambda i, j: self.scorePairwise(sequences[i], sequences[j])
        dist_vect = np.vectorize(dist_func)
        pre_matrix = np.fromfunction(dist_vect, shape = (len(sequences), len(sequences)))

        # scale data to between [0, 1]
        #min_max_fit = MinMaxScaler().fit_transform(pre_matrix)
        df_matrix = pd.DataFrame(pre_matrix, index = sequences, columns = sequences)
    
        return df_matrix

    def skbioMatrix(self):
        """
        test it out, yo.
        """

        pass


    def saveMatrix(self, save_as):
        """
        save distance matrix as specified file.
        """

        pass