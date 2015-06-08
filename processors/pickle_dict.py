__description__ = \
"""
Import a pickle_dict file.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-05-25"

import os, pickle, shutil
from . import BaseProcessor

class PickleDictProcessor(BaseProcessor):
    """
    Copy in a pickled dictionary file and load the thing into memory.
    """

    def process(self,pickle_file,data_label="raw_counts"):
        """
        Load in a pickled dictionary.
        """

        self.data_label = data_label

        internal_name = os.path.join(self.getProperty("expt_name"),
                                     "{:s}.pickle".format(self.data_label))

        file_label = "{:s}-file".format(self.data_label)
        self.addProperty(file_label,internal_name) 
        shutil.copy(pickle_file,self.getProperty(file_label))

        tmp_dict = pickle.load(open(self.getProperty(file_label),'rb'))

        self.addProperty(self.data_label,tmp_dict)

        self._do_not_write_to_json.append(self.data_label)

    @property
    def data(self):
        """
        Return the contents of the dictionary stored in the pickle dict.
        """

        return self.getProperty(self.data_label)
