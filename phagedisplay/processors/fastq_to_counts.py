__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-28"

import os, gzip, pickle, re

from . import BaseProcessor

class FastqSeqCounter:
    """
    Class for converting a set of fastq nucleotide sequences, collected over 
    multiple rounds of selection, and turning them into a single dictionary of the
    form:

    {"seq1":[10,20,100,...]}

    where the integers count the number of times that this sequence was seen in 
    each round.  
    """

    def __init__(self,bad_pattern="[*X]",phage_term="GGG*AET",seq_length=12):

        self._GENCODE = {
            'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
            'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
            'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
            'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
            'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
            'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
            'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
            'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
            'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
            'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
            'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
            'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
            'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
            'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
            'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
            'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

        self._bad_pattern = re.compile(bad_pattern)
        self._phage_term = phage_term
        self._seq_length = seq_length

    def _translate(self,sequence):
        """
        Translate a nucleotide sequence into a protein sequence.  If there is a
        problem, write an "X" into the sequence.
        """

        try:
            return "".join([self._GENCODE[sequence[3*i:3*i+3]]
                            for i in range(len(sequence)//3)])
        except KeyError:
            out = []
            for i in range(len(sequence)//3):
                try:
                    out.append(self._GENCODE[sequence[3*i:3*i+3]])
                except KeyError:
                    out.append("X")
            return "".join(out)

    def _qualityCheck(self,sequence):
        """ 
        Make sure that the C-terminus of the peptide is what we think it should
        be and that there aren't any "_bad_patterns" in the protein sequence. 
        (As of this writing, _bad_pattern will catch "X" and premature stop 
        codons). 
        """

        if sequence[self._seq_length:(self._seq_length+len(self._phage_term))] == self._phage_term:
            if not self._bad_pattern.search(sequence[:self._seq_length]):
                return True

        return False
   
    def _processSingleFile(self,fastq_file):
        """
        Create a set of good and bad pattern dicts for a given fastq file.
        """

        good_count_dict = {}
        bad_count_dict = {}

        get_line = False
        with gzip.open(fastq_file,'r+') as f:
       
            for l in f:
            
                l_ascii = l.decode("ascii")

                # If we see @, take the next line
                if l_ascii.startswith("@"):
                    get_line = True
                    continue

                # If get_line is False, continue.  Otherwise, set get_liune back
                # to False and move on to do the parsing business.
                if not get_line:
                    continue
                else:
                    get_line = False

                # Translate the sequence
                sequence = self._translate(l_ascii.strip())

                # Record it in either the good or bad dict, depending on its
                # quality score
                if self._qualityCheck(sequence):

                    key = sequence[0:self._seq_length]
                    try:
                        good_count_dict[key] += 1
                    except KeyError:
                        good_count_dict[key] = 1
                else:
                    try:
                        bad_count_dict[sequence] += 1
                    except KeyError:
                        bad_count_dict[sequence] = 1

        return good_count_dict, bad_count_dict

    def _compressDictSet(self,list_of_dicts):
        """
        Take a list of dictionaries, each corresponding to the counts for each
        peptide in a round, and create a single output dictionary keying 
        sequence to the number of counts in each round.  Any round that was 
        done, but not seen in the dicts, is given a "None" entry.
        """
    
        # all_keys has every sequence seen
        all_keys = []
        template = [None for i in range(len(list_of_dicts))]
        for i, a in enumerate(list_of_dicts):
            if a:
                template[i] = 0
                all_keys.extend(a.keys())

        # Create final dictionary that we'll populate with values below
        out_dict = dict([(a,template[:]) for a in all_keys])

        # Populate the final dictionary
        for i in range(len(list_of_dicts)):
            if not list_of_dicts[i]:
                continue

            for key, value in list_of_dicts[i].items():
                out_dict[key][i] = value

        return out_dict

    def processFastqFiles(self,fastq_file_list):
        """
        Take a set of fastq files and count them.
        """

        # Count good and bad reads and put them into lists of dictionaries
        all_good_dicts = [] 
        all_bad_dicts = [] 
        for f in fastq_file_list:

            good_counts = None
            bad_counts = None
            if f != None:
                print("Processing %s" % f)
                good_counts, bad_counts = self._processSingleFile(f)
            
            all_good_dicts.append(good_counts)
            all_bad_dicts.append(bad_counts)

        # Create a final dictionary for the good counts
        good_dict = self._compressDictSet(all_good_dicts)
    
        # Create a final dictionary for the bad counts
        bad_dict = self._compressDictSet(all_bad_dicts)

        return good_dict, bad_dict

    @property
    def bad_pattern(self):
        """
        Return regular expression we use to look for "badness"
        """
        return self._bad_pattern.pattern

    @property
    def phage_term(self):
        """
        Return the amino acid sequence of the C-terminus in the phage itself
        """
        return self._phage_term

    @property
    def seq_length(self):
        """
        Return expected length of peptide sequences.
        """
        return self._seq_length

class FastqToCountsProcessor(BaseProcessor):
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
        p = FastqSeqCounter()

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
        self.addProperty("good-counts",good_counts)
        self.addProperty("bad-counts",bad_counts)

        self._do_not_write_to_json.append("good-counts")
        self._do_not_write_to_json.append("bad-counts")

        self.saveFile(overwrite=True)


    @property
    def data(self):
        """
        """

        return self.getProperty("good-counts")
