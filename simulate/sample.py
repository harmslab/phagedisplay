# Stupidly simple sampling object... 
# This is a just an object wrapping numpy's sampling method `choice`
# The point is that we add properties and methods to SamplePool object that help
# us query the pool in simulations. 

import numpy as np

class SamplePool(object):
    
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
        
        
