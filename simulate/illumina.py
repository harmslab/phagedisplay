__description__ = \
"""
Sample from a pool in Illumina fashion.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-28"

import numpy as np

class IlluminaRun(object):
    
    def __init__(self,num_reads=100e6):
        """
        """

        self._num_reads = 100e6

    
    @property
    def num_reads(self):
        """
        Number of Illumina reads.
        """

        return self._num_reads
