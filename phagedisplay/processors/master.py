__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-28"

import os, gzip,  pickle

import phagedisplay
from phagedisplay.processors import BaseProcessor

class MasterProcessor(BaseProcessor):
    """
    Processor that will hold an entire phage display experiment.
    """

    def process(self,input_subprocessor=None,**kwargs):
        """
        """

        # No suprocessor inputs specified.  So either do nothing (no 
        # subprocessor around), take in outside data via **kwargs, or 
        # use the last subprocessor.
        if input_subprocessor == None:

            # No subprocessors, nothing to do
            if len(self._subprocessors) == 0:
                return None

            # Has only one subprocessor, so all data must come from outside
            elif len(self._subprocessors) == 1:
                self._subprocessors[-1].process(**kwargs)

            # Has more than one subprocessorBase.  Take data from previous
            # processor plus the kwargs
            else:
                self._subprocessors[-1].process(self._subprocessors[-2].data,
                                                **kwargs)

        # work on the subprocessor specified in the self.process call
        else:
            self._subprocessors[-1].process(input_subprocessor.data,
                                            kwargs) 
        

        self.saveFile(overwrite=True)

    @property
    def data(self):
        """
        """

        return self._subprocessors


