"""
To save memory and avoid big-ole dicts, treat sequences as base-20 numbers that
are converted to integers.  


"""

import numpy as np

class SamplePool(object):
    """ 
    Stupidly simple sampling object... 
    This is a just an object wrapping numpy's sampling method `choice`
    The point is that we add properties and methods to SamplePool object that
    help us query the pool in simulations. 
    """   
 
    def __init__(self, sequences, weights):
        """ """
        self._sequences = np.array(sequences)
        self._weights = np.array(weights)
        
    @property
    def sequences(self):
        """ Get sequences. """
        return self._sequences
        
    @property
    def weights(self):
        """ Get weights. """
        return self._weights
        
    def sample(self, size):
        """ Return a sample. """
        return np.random.choice(self._sequences, size=size, replace=True, p=self._weights)

    def amplify(self,observation):
        """
        Amplify samplified sequences into a set of sequences with associated 
        probability distributions.  As implemented now, this won't distort
        any frequencies, but this method can be replaced in a subclass.
        """

        # What we *hope* is happening in the lab...
        # amp_factor = 2**num_rounds
        # sequences = [s*amp_factor for s in sequences]
        
        counts = np.bincount(observation)
        values = np.nozero(counts)[0]
        weights = (1.0*counts[values])/sum(counts[values])
       
        return values, weights
 
        
class Pool:
    """
    Basic class for generating a pool from which to sample.  
    """ 

    def __init__(self,sequence_length=12,alphabet=(0,1)):
        """
        """

        self._sequence_length = sequence_length
        self._alphabet = alphabet[:]
  
    @property
    def sequence_length(self):
        """ Get sequence length. """
        return self._sequence_length

    @property
    def alphabet(self):
        """ Get alphabet. """
        return self._alphabet


    
