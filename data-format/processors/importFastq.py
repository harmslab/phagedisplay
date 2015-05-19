__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-28"

import os, gzip

class ImportFastq(base.ProcessorBase):
    """
    Class to import and hold onto a set of fastq files corresponding to a set
    of enrichment rounds.
    """

    def addFastqFile(self,file_name,file_round):
        """
        Add a fastq file to the experiment.  
        
        Args:
            file_name: name of data file
            file_round: experimental round data file corresponds to (int >= 0)
        """

        key = "fastq-file_{0:8d}".format(file_round)

        base_dir = self.getProperty("base_dir")
        value = os.path.join(base_dir,"{}.gz".format(os.path.split(file_name)[-1]))
       
        # Bring in the fastq file, compressing along the way 
        with open(file_name,'rb') as fastq_in:
            with gzip.open(value, 'wb') as gzipped_fastq_out:
                gzipped_fastq_out.writelines(fastq_in) 

        # Add the property keying to the newly created file
        self._addProp(key,value,blob.FileBlob) 


    def process(self,file_list,round_list=None):
        """
        Given a set of fastq files and round numbers, load them in.
        """   
 
        if not round_list:
            round_list = range(len(file_list))

        for i, f in enumerate(file_list):
            self.addFastqFile(f,round_list[i])

    def load(self):
        """
        Overwrite generic loader, adding extra sanity check to make sure that
        specified fastq files actually exist. 
        """
        
        super().load()

        # Make sure the loaded files actually exist.
        for k in self._props:
            if k.startswith("fastq-file_"):
                if os.path.exists(self.getProperty(k)):
                    err = "Could not find file: {:s}\n".format(self.getProperty(k))
                    raise ValueError(err)
        
    @property
    def data(self):
        """
        Return list of loaded fastq files, with "None" for rounds not done.
        """
  
        round_list = []
        file_list = [] 
        for k in self._props:
            if k.startswith("fastq-file_"):
                round_number = int(k.split("_")[1]
                filename = self.getProperty(k) 
             
                round_list.append(round_number)
                file_list.append(filename)

        out_list = [None for i in range(max(round_list))]
        for i range(len(round_list)):
            out_list[round_list[i]] = file_list[i]

        return out_list
