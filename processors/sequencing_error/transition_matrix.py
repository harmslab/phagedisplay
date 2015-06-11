__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-06-09"

import numpy as np

class TransitionMatrix:
    """
    """

    def __init__(self,base_value=1.0):
        
        self._alphabet = []
        self._transition_matrix = []
        self._alphabet_length = len(self.alphabet)

    def setup(self,base_value)

        self._base_value = base_value
        self._letter_to_index_dict = dict((a,i) for i, a in enumerate(self._alphabet))
        self._alphabet_length = len(self.alphabet) 
     
        self._transition_matrix = np.ones((self.alphabet_length,
                                           self.alphabet_length),dtype=float)

        # Set A->A = 0, T->T = 0, etc.
        self_to_self = np.array(zip(range(self.alphabet_length),range(self.alphabet_length)))
        self._transition_matrix[self_to_self] = 0.

        self._transition_matrix = self._transition_matrix*self._base_value

    @property
    def alphabet(self):
        return self._alphabet

    @property
    def transition_matrix(self):
        return self._transition_matrix

    @property
    def alphabet_length(self):
        return self._alphabet_length

    def letter_to_index(self,letter):
        """
        """

        return self._letter_to_index_dict[letter]

class FlatDNATM(TransitionMatrix):
    """
    """
    
    def __init__(self,base_value=1.0):
        
        self._alphabet = ["A","T","G","C"]
        self.setup(base_value)

class FlatProteinTM(TransitionMatrix):
    """
    """
    
    def __init__(self,base_value=1.0):
      
        warn = "This transition matrix models all amino acid subsitutions as "
        warn += "equally probable!  Really only for testing..."
        Warning(warn)
 
        self._alphabet = ["A","C","D","E","F","G","H","I","K","L","M","N","P",
                          "Q","R","S","T","V","W","Y"]
        self.setup(base_value)
