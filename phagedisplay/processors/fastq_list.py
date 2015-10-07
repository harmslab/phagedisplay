__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-01-09"

import re, gzip, os, shutil

from . import BaseProcessor

class FastqListProcessor(BaseProcessor):
    """
    Class to import and hold onto a set of fastq files corresponding to a set
    of enrichment rounds.
    """

    def process(self,file_list):
        """
        Given a set of fastq files and round numbers, load them in.
        """
   
        initial_fastq = [None for i in range(len(file_list))]
        self.addProperty("fastq-files",initial_fastq)
    
        expt_name = self.getProperty("expt_name")
        for i, f in enumerate(file_list):
            if f:
            
                self._logger("Importing {:s}".format(f))

                # If it's already compressed, just copy it in.
                if f[-3:] == ".gz":
                    value = os.path.join(expt_name,os.path.split(f)[-1])
                    shutil.copy(f,value)

                # Otherwise, compress it ourselves
                else:
                    value = os.path.join(expt_name,"{}.gz".format(os.path.split(f)[-1]))
 
                    # Bring in the fastq file, compressing along the way 
                    with open(f,'rb') as fastq_in:
                        with gzip.open(value, 'wb') as gzipped_fastq_out:
                            gzipped_fastq_out.writelines(fastq_in) 


                # Update the "None" entry in the fastq-files list to be the filename for
                # for this round.
                self.getProperty("fastq-files")[i] = value
        
        self.saveFile(overwrite=True)


    @property
    def data(self):
        """
        Return list of loaded fastq files, with "None" for rounds not done.
        """
 
        return self.getProperty("fastq-files")
