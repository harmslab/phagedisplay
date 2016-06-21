import pickle
import numpy as np

def get_logK(prev_counts,this_counts):
    """
    """
    
    if len(prev_counts) != len(this_counts):
        err = "Count arrays must have the same number of entries.\n"
        raise ValueError(err)
        
    num_prev_zero = np.sum(prev_counts == 0)
    num_this_zero = np.sum(this_counts == 0)
    if num_prev_zero != 0 or num_this_zero != 0:
        err = "All array entries must be non-zero\n"
        raise ValueError(err)
        
    # output array [[log(K), lower_bound,upper_bound],...]
    out = np.zeros((len(prev_counts),3),dtype=float)
    
    # create index matrix and its inverse
    # [[-1,0,0],
    #  [-1,1,0],
    #  [-1,0,1]]
    # column 0 is omega (subtracted from everyone).  
    X = np.eye(len(prev_counts))
    X[:,0] = -1 
    inv_X = np.linalg.inv(X)
    
    # Determine the enrichment ratio for each set of counts
    enrich_ratio = (this_counts/np.sum(this_counts)) * (np.sum(prev_counts)/prev_counts)

    # Use linear model to solve for log(K)
    out[:,0]= np.dot(inv_X,np.log(enrich_ratio))
    
    # Now propagate uncertainty
    sigma_ratio = enrich_ratio * np.sqrt(this_counts**(-3) + prev_counts**(-3))
    K_sigma = np.sqrt(np.dot(np.abs(inv_X),sigma_ratio**2))
   
    K = np.exp(out[:,0])
    out[:,1] = np.log(K - K_sigma)
    out[:,2] = np.log(K + K_sigma)
    
    # Grab value for omega
    omega = out[0,:]
    
    # Set log(K0) to zero.  This is our reference state
    out[0,0] = 0.0
    out[0,1] = 0.0
    out[0,2] = 0.0
    
    return out, omega

def pre_treat_arrays(seq_array,prev_counts,this_counts):
    """
    Get count arrays ready.
    """
    
    # Arrays of zero counts
    prev_zero = (prev_counts == 0)
    this_zero = (this_counts == 0)
    
    # Toss data for which both previous and current counts are 0    
    at_least_one_not_zero = np.logical_not(prev_zero * this_zero)
    
    seq_array   =   seq_array[at_least_one_not_zero]
    prev_counts = prev_counts[at_least_one_not_zero]
    this_counts = this_counts[at_least_one_not_zero]
    
    # Update zero count arrays to reflect toss
    prev_zero = prev_zero[at_least_one_not_zero]
    this_zero = this_zero[at_least_one_not_zero]
    
    # Collapse everyone for which previous or current is zero into a single large bin
    either_is_zero = prev_zero + this_zero
    
    prev_zero_counts = np.sum(prev_counts[either_is_zero])
    this_zero_counts = np.sum(this_counts[either_is_zero])
    
    # Now, remove entries for wich either is zero
    neither_is_zero = np.logical_not(either_is_zero)
    
    seq_array   =   seq_array[neither_is_zero]
    prev_counts = prev_counts[neither_is_zero]
    this_counts = this_counts[neither_is_zero]
    
    # Sort array from smallest to largest number of counts
    array_order = np.argsort(prev_counts + this_counts)
    
    seq_array   =   seq_array[array_order]
    prev_counts = prev_counts[array_order]
    this_counts = this_counts[array_order]
    
    # Add final junk category to the end of the array
    if prev_zero_counts > 0 and this_zero_counts > 0:
        seq_array =   np.append(seq_array,"".join(["X" for i in range(len(seq_array[0]))]))
        prev_counts = np.append(prev_counts,prev_zero_counts)
        this_counts = np.append(this_counts,this_zero_counts)
    
    return seq_array, prev_counts, this_counts

class Stuff:

    def __init__(self,delta_file,pickled_dict):
        """
        """
        
        self._delta_file = delta_file
        self._pickled_dict = pickled_dict
        
        self._load_delta_file()
        self._load_count_dict()
        
        self._filter = None

    
    def _load_delta_file(self):

        self._seq = []
        self._reference_df = []
        self._competitor_df = []
        
        with open(self._delta_file) as f:
            for line in f:
                col = line.split()

                if col[0] == "seq" or col[1] == "nan" or col[2] == "nan":
                    continue

                try:
                    self._seq.append(col[0].strip())
                    self._reference_df.append(float(col[1]))
                    self._competitor_df.append(float(col[2]))
                except ValueError:
                    continue
                
        
        self._seq_dict = dict([(s,i) for i,s in enumerate(self._seq)])
        
        self._seq = np.array(self._seq)
        self._reference_df = np.array(self._reference_df)
        self._competitor_df = np.array(self._competitor_df)
        

    def _load_count_dict(self):
        """
        """
        
        self._count_dict = pickle.load(open(self._pickled_dict,"rb"))
        
        self._count_seq = list(self._count_dict.keys())
        self._count_seq_index = dict([(k,i) for i, k in enumerate(self._count_seq)])
        
        self._num_count_rounds = len(self._count_dict[list(self._count_dict.keys())[0]])
        self._num_count_seq = len(self._count_seq)
        
        self._count_array = np.zeros((self._num_count_seq,self._num_count_rounds),dtype=int)
        
        for i, k in enumerate(self._count_seq):
            self._count_array[i,:] = np.array(self._count_dict[k])
            
        
    def apply_filter(self,ref_df_cutoff=None,comp_df_cutoff=None):
        """
        Filter data by some criteria.  Filters are non-destructive and can be 
        removed by remove_filter method.
            
            ref_df_cutoff: take sequences in which change in frequency for the
                           reference change is > ref_df_cutoff.
            comp_df_cutoff: take sequences in which change in frequency for the
                            competitor peptide is < comp_df_cutoff.
         
        """
        
        if ref_df_cutoff != None: 
            meets_ref_filter = self._reference_df > ref_df_cutoff
        else:
            meets_ref_filter = 1

        if comp_df_cutoff != None:
            meets_comp_filter = self._competitor_df < comp_df_cutoff
        else:
            meets_comp_filter = 1
        
        self._filter = meets_ref_filter * meets_comp_filter
        
    def remove_filter(self):
        """
        Wipe out the filter.
        """
        
        self._filter = None

    def round_counts(self,round_number=1):
        """
        Return the counts ofr a given round.
        """
        
        if round_number >= self._num_count_rounds:
            err = "round not found\n"
            raise ValueError(err)
                    
        indexes = np.array([self._count_seq_index[s] for s in self.sequences])
        
        return self._count_array[indexes,round_number]

        
    @property
    def sequences(self):
        
        if self._filter != None:
            return self._seq[self._filter]
        else:
            return self._seq
       
X = Stuff("p000-p100_r1-r2_peptide-delta.txt",
          "processed-experiments/NCX1_0_CaE/raw-counts/good-counts.pickle")

seq_array = X.sequences
prev_round = X.round_counts(1)
this_round = X.round_counts(2)

new_seq_array, new_prev, new_this = pre_treat_arrays(seq_array,prev_round,this_round)
logK, omega = get_logK(new_prev,new_this)
