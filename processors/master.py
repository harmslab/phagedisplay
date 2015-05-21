__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-28"

import os, gzip,  pickle
from . import BaseProcessor

class MasterProcessor(BaseProcessor):
    """
    Processor that will hold an entire phage display experiment.
    """

    def process(self,**kwargs):
        """
        """

        # No subprocessors, nothing to do
        if len(self._subprocessors) == 0:
            return None

        # Has only one subprocessor, so all data must come from outside
        elif len(self._subprocessors) == 1:
            self._subprocessors[-1].process(**kwargs)

        # Has more than one subprocessorBase.  Take data from previous processor,
        # plust the kwargs
        else:
            self._subprocessors[-1].process(self._subprocessors[-2].data,**kwargs)

        self.saveFile()

    @property
    def data(self):
        """
        """

        return self._subprocessors


