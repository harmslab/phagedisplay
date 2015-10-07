__description__ = \
"""
Use a maximum likelihood model to estimate the actual frequencies of clones in a
library given that there are sequencing errors.  Follows Lynch paper.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-06-09"

import numpy as np
#from . import BaseProcessor
import itertools

# -----------------------------------------------------------------------------

def find_seq_neighbors(seq,max_num_mutations=2,alphabet=("A","T","G","C")):
    """
    Take a sequence and generate a list of all possible sequence neighbors 
    within max_num_mutations.  

        input:

        seq: sequence string (assumes all letters are within alphabet)
        max_num_mutations: integer indicating how many mutations away to walk
        alphabet: possible states at each site.  

        output:
    
        all_possible_neighbors: a list of strings containing all possible 
                                neighbors to seq.
    """
        
    wt_seq = list(seq)
    num_sites = len(wt_seq)
    all_possible_neighbors = []

    # For all possible numbers of mutations (0 through max_num...)
    for i in range(max_num_mutations+1):

        # Create a list of all possible combinations of i states given the
        # alphabet
        state_combos = list(itertools.combinations_with_replacement(alphabet,(i)))

        # Go through all possible "num_sites choose i" combinatons of site indexes
        for sites in itertools.combinations(range(num_sites),i):

            # Now go through possible state combinations for these sites.
            for j in range(len(state_combos)):

                # Do the actual mutations to the string
                mutated_seq = wt_seq[:]
                for k in range(i):
                    mutated_seq[sites[k]] = state_combos[j][k]

                all_possible_neighbors.append("".join(mutated_seq))

    return all_possible_neighbors
       


class Rocker:

    def __init__(self,count_dict,num_steps,alphabet):
        """
        """

        self.num_steps = num_steps
        self.alphabet = alphabet

        seq_list = count_dict.keys()

        zero_tuple = tuple(0 for c in count_dict[seq_list[0]] if c != None)
        seq_dict = {}
        for i in range(len(seq_list)):
            seq_dict[seq_as_string[i]] = tuple([c for c in count_dict[seq_list[i]] if c != None])

        neighbors = np.zeros((num_seq,num_neighbors,num_rounds),dtype=float)
        for i, s in enumerate(sequences):
            for j, neighbor in enumerate(find_seq_neighbors(s,num_steps,self.alphabet)):

                # If we've seen the sequence and have counts for it, record. 
                if seq_dict[s2]:
                    neighors[i,j,:] = seq_dict[s2]

                # Otherwise, set counts to 0.
                else:
                    neighbors[i,j,:] = zero_tuple
        
        # Guesses
        self.param_guess = np.zeros((self.num_patterns*2),dtype=float) 
        for i in range(self.num_patterns):
            y = self.y_obs[i,:]
            conc_guess = self.y_obs[i,0]
            K_guess = 1.0 

            try:
                indep_param, cov = curve_fit(self._indepFit,self.rounds,y,
                                             p0=(conc_guess,K_guess),maxfev=10000)

                if indep_param[0] <= 0 or indep_param[1] <= 0:
                    raise RuntimeError
                self.param_guess[i] = np.log(indep_param[0]) 
                self.param_guess[self.num_patterns + i] = np.log(indep_param[1])
            except RuntimeError:
                self.param_guess[i] = np.log(1/self.num_patterns) 
                self.param_guess[self.num_patterns + i] = np.log(1.0)

        # Create bound list.  conc must be between 0 and 1, K must be positive
        self.bounds = [(None,None) for i in range(self.num_patterns)]   
        self.bounds.extend([(None,None) for i in range(self.num_patterns)])

        self.bounds[self.num_patterns] = (-1.0,-1.0)
        self.param_guess[self.num_patterns] = -1.0

    def runRegression(self,maxiter=100000000000):
        """
        Run a regression of the objective function given our bound, constraints,
        etc.
        """

        self.fit_result = minimize(fun=self._objective,x0=self.param_guess,
                                   bounds=self.bounds,
                                   options={"maxiter":maxiter,"maxfun":maxiter})

    def returnParam(self):
        """
        Return the fit parameters.  
        """

        initial_theta = self.fit_result.x[:self.num_patterns]
        K_values = self.fit_result.x[self.num_patterns:]

        out = np.zeros((len(K_values),3),dtype=float)
        for i in range(3):
            out[:,i] = np.exp(self.log_degeneracy + initial_theta + K_values*(i+1))
            out[:,i] = out[:,i]/sum(out[:,i])

        return initial_theta, K_values, out   

    def _objective(self,param):
        """
        Objective function to minimize.  Goes like:
        """

        for j, n in enumerate(self.rounds):
            self.p = np.exp(self.log_degeneracy + param[:self.num_patterns] + param[self.num_patterns:]*n)
            self.y_calc[:,j] = self.p/np.sum(self.p)

        return np.sum(np.power(self.y_calc-self.y_obs,2))

"""
class SequencingErrorProcessor(BaseProcessor):

    def process(self):

        pass

    @property
    def data(self):

        pass

"""
