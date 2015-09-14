__description__ = \
"""
Use a maximum likelihood model to estimate the actual frequencies of clones in a
library given that there are sequencing errors.  Follows Lynch paper.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-06-09"

import numpy as np
#from . import BaseProcessor

def _possibleNeighborsRecursive(sequence,num_steps,mutations,alphabet):
    """
    Generate a list of unique neighbors within num_steps hamming distance of the
    starting sequence given some alphabet.
    """
    
    neighbors = []
    original_sequence = list(sequence)
    original_mutations = mutations[:]

    # Go through each site in the sequence
    for i in range(len(original_sequence)):

        # At each site, go through each letter in the alphabet
        for a in alphabet:

            # Make the change and append to neighbors
            new_sequence = original_sequence[:]
            if new_sequence[i] == a:
                continue

            new_sequence[i] = a
            joined_sequence = "".join(new_sequence)
       
            mutations = original_mutations[:]
            mutations.append((original_sequence[i],i,a))

            neighbors.append((joined_sequence,tuple(mutations)))
           
            # Recursive call to take care of neighbors of neighbors 
            if num_steps > 1:
                neighbors.extend(_possibleNeighborsRecursive(joined_sequence,num_steps-1,mutations,alphabet))


    # Take unique sequences
    neighbors = list(dict([(n,[]) for n in neighbors]).keys())
   
    return neighbors 
             
def possibleNeighbors(sequence,num_steps,alphabet):
    """
    Generate a list of unique neighbors within num_steps hamming distance of the
    starting sequence given some alphabet.
    """
    
    neighbors = _possibleNeighborsRecursive(sequence,num_steps,[],alphabet)

    final_neighbors = {}
    for i in range(len(neighbors)):

        neighbor_sequence = neighbors[i][0]
        mutations = neighbors[i][1]

        if neighbor_sequence == sequence:
            continue

        try:
            current = final_neighbors[neighbor_sequence]
            if len(current) > len(mutations):
                final_neighbors[neighbor_sequence] = mutations[:]
        except KeyError:
            final_neighbors[neighbor_sequence] = mutations[:]


    for f in final_neighbors:
        print(f,final_neighbors[f],len(final_neighbors[f]))
    print(len(final_neighbors))

    return final_neighbors


possibleNeighbors("YYYYYYYYYYYY",2,["A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y"])

class Rocker:

    def __init__(self,count_dict,num_steps,alphabet):
        """
        """

        self.num_steps = num_steps
        self.alphabet = alphabet

        seq_list = count_dict.keys()
        seq_as_string = [converter(s) for s in seq_list]

        #round_list = [SOMEHOW POPULATE HERE]

        seq_dict = {}
        for i in range(len(seq_list)):
            seq_dict[seq_as_string[i]] = tuple([c for c in count_dict[seq_list[i]] if c != None])

        neighbors = np.zeros((num_seq,num_neighbors,num_rounds),dtype=float)
        for i, s in enumerate(sequences):
            for j, neighbor in enumerate(possibleNeighbors(s,num_steps)):

                # If we've seen the sequence and have counts for it, record. 
                if seq_dict[s2]:
                    neighors[i,j,:] = seq_dict[s2]

                # Otherwise, set counts to 0.
                else:
                    neighbors[i,j,:] = 0
        

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
