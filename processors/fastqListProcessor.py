__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-01-09"

import re, gzip, os
from processors import processorBase

class FastqListProcessor(processorBase.ProcessorParent):
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

                value = os.path.join(expt_name,"{}.gz".format(os.path.split(f)[-1]))
      
                print("Importing and compressing {:s}".format(f))
 
                # Bring in the fastq file, compressing along the way 
                with open(f,'rb') as fastq_in:
                    with gzip.open(value, 'wb') as gzipped_fastq_out:
                        gzipped_fastq_out.writelines(fastq_in) 

                # Update the "None" entry in the fastq-files list to be the filename for
                # for this round.
                self.getProperty("fastq-files")[i] = value

    @property
    def data(self):
        """
        Return list of loaded fastq files, with "None" for rounds not done.
        """
 
        return self.getProperty("fastq-files")
