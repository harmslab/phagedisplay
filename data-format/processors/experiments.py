__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-28"

import os, gzip, container, pickle

import processFastq

class TotalContainer(container.ContainerParent):
    """
    Container that will hold an entire phage display experiment.
    """
    
    def process(self,**kwargs):
        """
        """

        # No subcontainers, nothing to do
        if len(self._subcontainers) == 0:
            return None

        # Has only one subcontainer, so all data must come from outside
        elif len(self._subcontainers) == 1:
            self._subcontainers[-1].process(**kwargs)

        # Has more than one subcontainer.  Take data from previous container,
        # plust the kwargs
        else:
            self._subcontainers[-1].process(self._subcontainers[-2].data,**kwargs)

        self.saveFile()


class FastqListContainer(container.ContainerParent):
    """
    Class to import and hold onto a set of fastq files corresponding to a set
    of enrichment rounds.
    """

    def process(self,file_list):
        """
        Given a set of fastq files and round numbers, load them in.
        """
   
        initial_fastq = [None for i in range(len(file_list))]
        self._addProperty("fastq-files",initial_fastq)
    
        for i, f in enumerate(file_list):
            if f:
                self._addFastqFile(f,i)

    def _addFastqFile(self,file_name,file_round):
        """
        Add a fastq file to the experiment.  
        
        Args:
            file_name: name of data file
            file_round: round of experiment in file
        """

        expt_name = self.getProperty("expt_name")
        value = os.path.join(expt_name,"{}.gz".format(os.path.split(file_name)[-1]))
      
        print("Importing and compressing {:s}".format(file_name))
 
        # Bring in the fastq file, compressing along the way 
        with open(file_name,'rb') as fastq_in:
            with gzip.open(value, 'wb') as gzipped_fastq_out:
                gzipped_fastq_out.writelines(fastq_in) 

        # Update the "None" entry in the fastq-files list to be the filename for
        # for this round.
        self.getProperty("fastq-files")[file_round] = value

    @property
    def data(self):
        """
        Return list of loaded fastq files, with "None" for rounds not done.
        """
 
        return self.getProperty("fastq-files")

class FastqToCountsContainer(container.ContainerParent):
    """
    """

    def process(self,fastq_filenames,
                good_counts_pickle=None,
                bad_counts_pickle=None):
        """
        """

        if not isinstance(fastq_filenames,list):
            err = "process requires list of fastq filenames\n"
            raise ValueError(err)

        # Create object to count the fastq sequences
        p = processFastq.FastqSeqCounter()

        # Count
        good_counts, bad_counts = p.processFastqFiles(fastq_filenames)

        # Write out pickle files 
        if not good_counts_pickle:
            good_counts_pickle = os.path.join(self.getProperty("expt_name"),
                                              "good-counts.pickle") 
        if not bad_counts_pickle:
            bad_counts_pickle = os.path.join(self.getProperty("expt_name"),
                                             "bad-counts.pickle") 
        f = open(good_counts_pickle,'wb')
        pickle.dump(good_counts,f)
        f.close()
     
        f = open(bad_counts_pickle,'wb')
        pickle.dump(bad_counts,f)
        f.close()

        # Record that we have the counts
        self._addProperty("good-counts",good_counts)
        self._addProperty("bad-counts",bad_counts)

        self._do_not_write_to_json.append("good-counts")
        self._do_not_write_to_json.append("bad-counts")

    @property
    def data(self):
        """
        """

        return self.getProperty("good-counts")
