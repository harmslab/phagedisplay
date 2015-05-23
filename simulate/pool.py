__description__ = \
"""
Classes for storing sequence pools that can then be sampled.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-25"

import string
import numpy as np

import utility

AMINO_ACIDS = ("A","C","D","E","F",
               "G","H","I","K","L",
               "M","N","P","Q","R",
               "S","T","V","W","Y")

BASES = ("A","T","G","C")

BINARY = (0,1)

class SeqIntegerMapper:
    """
    To save memory and avoid big-ole dicts, treat sequences as base-"alphabet
    size" numbers that are converted to base 10 integers.  Can handle alphabets
    up to 36 letters long. As written, we'll run out of long ints for 64 bit
    systems for 15+ amino acid peptides.  
    """
    
    def __init__(self,seq_length,alphabet=AMINO_ACIDS):
        """
        Initialize instance of class.
        """
        
        self._alphabet = ["{:s}".format(s) for s in alphabet]
        self._base = len(alphabet)
        self._seq_length = seq_length
        self.possible_digits = string.digits + string.ascii_lowercase
        
        self.max_int = self._base**(self._seq_length)
        
        self.string_to_base = dict([(letter,self.possible_digits[i])
                                    for i, letter in enumerate(alphabet)])
                
    @property
    def alphabet(self):
        """
        Get alphabet.
        """
        return self._alphabet

    @property
    def base(self):
        """
        Get base.
        """
        return self._base
    
    @property
    def seq_length(self):
        """ Get base. """
        return self._seq_length
        
    def seqToInt(self,sequence_string):
        """
        Return the base 10 integer equivalent of a sequence.
        """
    
        sequence_in_base = "".join([self.string_to_base[s]
                                    for s in sequence_string])
        
        return int(sequence_in_base,self._base)
    
    
    def intToSeq(self,sequence_integer):
        """
        Return the sequence encoded by this base 10 integer.
        """
            
        digits = list(self._alphabet[0])*self._seq_length
        i = 0
        while sequence_integer:
            digits[i] = self._alphabet[sequence_integer % self._base]
            sequence_integer //= self._base
            i += 1
          
        digits.reverse()
        
        return "".join(digits)

class Pool:
    """
    Class to hold a pool of sequences. 
    
    Important attributes:
   
         LISTS WHERE EACH ENTRY IS A SEQUENCE 
         _all_seq holds integer representations of every sequence in the initial
                        pool.
         _affinities    holds the relative affinities of each sequence in
                        _all_seq

         LISTS WHERE EACH ENTRY IS A ROUND
         _contents      is a list of arrays of integers ranging from 0 to
                        len(_all_seq)-1. Each array in the list corresponds to
                        one round in the experiment. Each integer in the arrays
                        maps back to the sequences stored in _all_seq (and the
                        affinities in _affinities)
         _counts        holds a list of integer arrays storing the counts of
                        each sequence in the pool for a given round. 
         _checkpoints   is a list of boolean variables that records whether this
                        particular round is interesting.  This is set using
                        addNewStep and is useful for keeping track of steps that
                        correspond to, say, experimental outputs.

    These private attributes can be accessed by a variety public functions. 
    
    The class also contains methods for generating pools and adding new rounds
    (updating the _contents and _counts lists).  
    """ 

    def __init__(self,sequence_length=12,alphabet=BINARY):
        """
        Initialize the class with basic pool properties.
        
        Args: (sequence_length) length of sequences in pool
              (alphabet) possible states at each site in the sequence
        """

        self.sequence_length = sequence_length
        self.alphabet = alphabet[:]
        
        # mapper allows us to convert between internal integer representation and
        # human-readable sequence representations.
        self.mapper = SeqIntegerMapper(self.sequence_length,
                                       tuple(self.alphabet[:]))
        
        # At this point, there are no sequences in the pool...
        self._pool_exists = False
    

    def createUniformPool(self,initial_pool_size,max_K=1e6):
        """
        Generate a pool of "initial_pool_size" sequences, sampling affinities
        randomly from 1 to max_K. Affinities are assigned using randomly chosen
        log(K) values, but affinity is stored as K values. 
        """
        
        # Create random initial sequences and count them.  mapper.max_int is the 
        # highest possible sequence.  For example, for a 5-mer peptide, this will
        # generate integer equivalents ranging from AAAAA to YYYYY.  
        initial_sample = np.random.randint(0,self.mapper.max_int,
                                           initial_pool_size)  
        content, counts = utility.uniqueCounter(initial_sample)
        
        # All sequences seen, as well as their srelative affinities
        self._all_seq = content
        self._affinities = 10**(np.random.uniform(low=0.0,high=np.log10(max_K),
                                                  size=self._all_seq.size))
        
        # Original pool of sequences
        self._contents = [np.array(range(self._all_seq.size))]
        self._counts = [counts]
        self._checkpoints = [True]
       
        # Now we have a pool to work with.  
        self._pool_exists = True
        
    def createSkewedPool(self,initial_pool_size,max_K=1e6,scale_reps=3):
        """
        """
        
        # Create random initial sequences and count them.  mapper.max_int is the 
        # highest possible sequence.  For example, for a 5-mer peptide, this will
        # generate integer equivalents ranging from AAAAA to YYYYY.  
        initial_sample = np.random.randint(0,self.mapper.max_int,
                                           initial_pool_size)  
        content, counts = utility.uniqueCounter(initial_sample)
        
        # All sequences seen, as well as their srelative affinities
        self._all_seq = content

        s = 1
        for i in range(scale_reps):
            s = s*np.random.uniform(low=0.0,high=1.,size=self._all_seq.size)
            #b = np.random.uniform(low=0.0,high=1.,size=self._all_seq.size)
            #c = 1   #np.random.uniform(low=0.0,high=1.,size=self._all_seq.size)
        #self._affinities = 10**(a*b*c*np.log10(max_K))
        self._affinities = 10**(s*np.log10(max_K))
        
        # Original pool of sequences
        self._contents = [np.array(range(self._all_seq.size))]
        self._counts = [counts]
        self._checkpoints = [True]
       
        # Now we have a pool to work with.  
        self._pool_exists = True
    
    def addNewStep(self,new_contents,new_counts,checkpoint):
        """
        Append a new round of selection.  

        FUTURE NOTE: Could make this cleverly write out current contents as 
        a pickle and then replace to allow large simulations. 
        """

        self._contents.append(new_contents)
        self._counts.append(new_counts)
        self._checkpoints.append(checkpoint)
   
    def reset(self):
        """
        Reset the pool to its initial state.
        """

        self._contents = self._contents[:1]
        self._counts = self._counts[:1]
        self._checkpoints = self._checkpoints[:1]
 
    def prettyPrint(self):
        """
        Print out the current round in tabular fashion.  
        """
        
        freq = self.current_counts/sum(self.current_counts)
        for i in self.current_contents:
            print(self.mapper.intToSeq(self._all_seq[self.current_contents[i]]),
                  self.current_counts[i],
                  self.current_affinities[i],
                  freq[i])  
    
    @property
    def pool_exists(self):
        """
        Has the pool actually been populated?
        """
        
        return self._pool_exists
    
    @property
    def current_contents(self):
        """
        The current contents of the pool. 
        """
        
        return self._contents[-1]
    
    @property
    def current_counts(self):
        """
        The current counts of all sequences in the current pool.  
        """
        return self._counts[-1]  

    @property
    def current_affinities(self):
        """
        The affinities of all of the sequences in the currrent pool. 
        """
        
        return self._affinities[self._contents[-1]]
   
    @property
    def checkpoints(self):
        """
        Checkpoint status of each round.
        """

        return self._checkpoints
 
    @property
    def all_seq(self):
        """
        All sequences ever seen in the pool.  
        """
        
        return self._all_seq
    
    @property
    def all_affinities(self):
        """
        The affinities of every sequence ever seen in the pool. 
        """
        
        return self._affinities 

    def round_contents(self,round_number):
        """
        Get contents of pool at round round_number.
        """
    
        return self._contents[round_number]
    
    def round_counts(self,round_number):
        """
        Get counts of sequences in pool at round round_number.
        """
    
        return self._counts[round_number]
    
    def round_affinities(self,round_number):
        """
        Get affinities of sequences in pool at round round_number.
        """
    
        return self._affinities[self._contents[round_number]]

    @property
    def round_count_dict(self):
        """
        Return a count dict in the same format that the experimental analysis 
        pipeline does:
            {"SEQUENCE":[round0_count,round1_count,round2_count...],...}
        """

        out_dict = {}
        rounds_to_write = [i for i, c in enumerate(self._checkpoints) if c]
        template = [0 for i in range(len(rounds_to_write))]
        for a in self._all_seq:
            seq = self.mapper.intToSeq(a)
            out_dict[seq] = template[:]

        for i, r in enumerate(rounds_to_write):
            contents = self.round_contents(r)
            counts = self.round_counts(r)

            for j in range(len(contents)):
                key = self.mapper.intToSeq(self._all_seq[contents[j]])
                out_dict[key][i] = counts[j]
            
        return out_dict 

       
