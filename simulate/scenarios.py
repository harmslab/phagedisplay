
import sys, pickle
from simulate import pool, samplers, utility, parameters
from scipy.optimize import minimize

class StandardExperiment:
    """
    A standard, multiple round phage display experiment simulation.  
    """
    
    def __init__(self,sequence_length=8,alphabet=pool.AMINO_ACIDS):
        """
        Initialize the run.
        """

        # Calculate number of possible sequences
        self.sequence_length = sequence_length
        self.alphabet = alphabet
        self.total_possible = len(alphabet)**sequence_length
        
        self._input_pool_dict = {}
        self._actual_out_dict = {}
        self._illumina_out_dict = {}

    def create(self,param_set=parameters.lwheeler00,affinity_max=1e6,skew=3,conc_const_guess=10.,specified_conc_const=-1):
        """
        Create a simulation engine with a concentration constant that reproduces the
        observed experimental value.

        affinity_max: affnity of the maximum sequence (all sequences will have affinities
                      between 1 and affinity_max, so this really sets the dynamic range 
                      in the library).
        skew:         how much to skew the library towards low affinity binders.  The 
                      higher the number, the greater the skew to low affinity.
        conc_const_guess: guess for starting concentration constant for optimization 
                      of this value.
        specified_conc_const: specify a value of the concentration constant directly.  
                      Negative values indicate that it should be ooptimized. 
        """
       
        for k in param_set.__dict__.keys():
            if not k.startswith("_"):
                setattr(self,k,getattr(param_set,k))
       
        try: 
            # Determine the number of sequences in each 
            self.company_tube = int(self.company_tube_fx*self.total_possible)
            self.pipet_onto_plate = int(self.pipet_onto_plate_fx*self.total_possible)
            self.protein_on_plate = int(self.protein_on_plate_fx*self.total_possible)
            self.amplify_post_bind = int(self.amplify_post_bind_fx*self.total_possible)
            self.illumina_reads = int(self.illumina_reads_fx*self.total_possible)
            self.first_round_recov = int(self.first_round_recov_fx*self.total_possible)
        except KeyError:
            err = "not all required parameters in specified param_dict\n"
            raise ValueError(err)

        # Create pipet sampler and amplifier.
        self.pipet = samplers.PipetteSampler()
        self.amplify = samplers.PhageAmplificationSampler()

        # Create a pool
        self.p = pool.Pool(sequence_length=self.sequence_length,alphabet=self.alphabet)
        self.p.createSkewedPool(self.pipet_onto_plate,affinity_max,scale_reps=skew)
        
        if specified_conc_const >= 0:
            self.conc_const = specified_conc_const
        else:
            # Optimize the concentration constant so, given the pool size, we reproduce what
            # we saw experimentally recovered on our first round.
            sys.stdout.write("# Optimizing concentration constant\n")
            conc_constant = minimize(self._conc_const_objective,conc_const_guess,
                                     method='nelder-mead',
                                     options={'xtol': 1, 'disp': False})
        
            sys.stdout.write("# Done\n")
            self.conc_const = conc_constant.x
        
        #Create a binding sampler
        self.screen = samplers.BindingSampler(conc_constant=10**self.conc_const)
        
        #Create an illumina sampler
        self.illumina = samplers.IlluminaSampler()
        
    def _conc_const_objective(self,conc_const):
        """
        Objective function for optimizing the concentration constant.
        """
                                                    
        screen = samplers.BindingSampler(conc_constant=(10**conc_const))
        self.pipet.runExperiment(self.p,self.pipet_onto_plate)
        recovery = screen.runExperiment(self.p,self.protein_on_plate,
                                        return_only_sample_size=True)
        self.p.reset()
       
        sys.stdout.write("# Current: {:.3f} {:d}. Target: {:d}\n".format(conc_const[0],
                                                                         recovery,
                                                                         self.first_round_recov))
        sys.stdout.flush()        
 
        return abs(self.first_round_recov - recovery)

    def run(self,num_rounds=3):
        """
        Run set of sampling rounds.
        """

        # Do sampling runs for num_rounds or until there are no phage left.
        try:
            print("{:10s}{:15s}{:15s}".format("round","num_unique","total_counts"))
            print("{:10d}{:15d}{:15d}".format(0,self.p.current_counts.size,sum(self.p.current_counts)))
            
            for i in range(num_rounds):
                self.pipet.runExperiment(self.p,self.pipet_onto_plate)
                self.screen.runExperiment(self.p,self.protein_on_plate)
                self.amplify.runExperiment(self.p,self.amplify_post_bind,checkpoint=True)
            
                print("{:10d}{:15d}{:15d}".format(i+1,self.p.current_counts.size,sum(self.p.current_counts)))
    
        except samplers.PoolIsEmptyError:
            pass

        # Input pool summary
        self._input_pool_dict = {}
        seqs = [self.p.mapper.intToSeq(s) for s in self.p.all_seq]
        affin = self.p.all_affinities
        initial_counts = self.p.round_counts(0)
        for i, s in enumerate(seqs):
            self._input_pool_dict[s] = (affin[i],initial_counts[i])
        
        # Result of sampling protocol, no illumina sampling
        self._actual_out_dict = self.p.round_count_dict
        
        # --------------------------------------------------------------------
        # Sample from the runs via simulated illumina sequencing of each round
        # --------------------------------------------------------------------
        
        # List of illumina samples to take
        to_take = [i for i,s in enumerate(self.p.checkpoints) if s]       
        
        # Sample the illumina counts at each round
        samples = []
        for t in to_take:
            contents, counts = self.illumina.runExperiment(self.p,
                                                           sample_size=self.illumina_reads,
                                                           round_to_sample=t)
            samples.append(list(zip(contents,counts)))
        
        # Create a list of unique sequences
        unique_seqs = []
        for s in samples:
            unique_seqs.extend([k[0] for k in s])
        unique_seqs = list(dict([(s,()) for s in unique_seqs]).keys())
        
        # Create a dictionary that will store sequence mapped to counts at each round
        template  = [0 for i in range(len(samples))]
        out_dict = {}
        for u in unique_seqs:
            out_dict[u] = template[:]
        
        # Append counts at each round for each sequence
        for i, s in enumerate(samples):
            for c in s:
                out_dict[c[0]][i] = c[1]

        # Convert the sequences to pretty, actual sequences
        for k in out_dict.keys():
            self._illumina_out_dict[self.p.mapper.intToSeq(self.p.all_seq[k])] = out_dict[k][:]
   
    def save(self,filename=None):
        """
        Write out the whole data structure to a pickle.
        """

        f = open(filename,'wb')
        pickle.dump(self.__dict__,f)
        f.close() 

    def load(self,filename):
        """
        Read a pickle file into an instance of the class.  (This will overwrite
        anything already stored in the instance that overlaps in __dict__, but
        keep the non-overlapping stuff).
        """

        f = open(filename,'rb')
        tmp_dict = pickle.load(f)
        f.close()

        self.__dict__.update(tmp_dict)

         
    @property
    def input_pool(self):
        """
        Give the input pool of sequences.  Sequences key to (affinity,initial_counts)
        """
        return self._input_pool_dict

    @property
    def actual_out(self):
        """
        Raw count output from rounds.
        """
        return self._actual_out_dict

    @property
    def illumina_out(self):
        """
        Count output, sampled via illumina.  
        """       
        return self._illumina_out_dict
    
    
