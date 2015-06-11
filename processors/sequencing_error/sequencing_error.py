__description__ = \
"""
Use a maximum likelihood model to estimate the actual frequencies of clones in a
library given that there are sequencing errors.  Follows Lynch paper.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-06-09"

import numpy as np
from scipy import misc
#from . import BaseProcessor

import transition_matrix

def _possibleNeighborsRecursive(sequence,num_steps,mutations,alphabet):
    """
    Recursive function to generate a list of neighbors.  Should actually call
    possibleNeighbors as this function will return duplicates of the sort
        CAA (one mutation A0C), CAA (two mutations, A0T, T0C).
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

    return final_neighbors.items()



class Rocker:

    def __init__(self,count_dict,num_steps,transition_matri):
        """
        """

        self.num_steps = num_steps
        self.tm = transition_matrix

        # Convert count_dictionary to a list, sorted from highest to lowest counts
        tmp_seq_list = list((v[1],v[0]) for v in count_dict.items())
        tmp_seq_list.sort(reverse=True)

        # Create core mapping arrays/dictionaries to link sequences back to their
        # indexes, counts, etc.
        self.seq_list =   [s[1] for s in tmp_seq_list]
        self.seq_counts = [s[0] for s in tmp_seq_list]
        self.seq_index_dict = dict((s,i) for i,s in enumerate(self.seq_list))
        self.num_real_seq = len(self.seq_list)

        # Add dummy sequences for each possible mutation combination
        possible_mutations = []
        for i, a in enumerate(self.tm.alphabet):
            for j, b in enumerate(self.tm.alphabet):
                if a != b:
                    possible_mutations.append((a,b))

        all_muts = []
        for i in range(num_steps):
            all_muts.extend(itertools.combinations(possible_mutations,i+1))

        self.seq_list.extend(all_muts)
        self.seq_counts.extend([0 for i in range(len(all_muts))])
        self.seq_index_dict.update((a,i+self.num_real_seq) for i, a in enumerate(all_muts))


        # Figure out how many neighbors are going to be possible for each sequence
        seq_length = len(self.seq_list[0])
        num_neighbors = 0
        for i in range(len(self.num_steps)):
            num_neighbors += misc.comb(seq_length,i)*(self.tm.alphabet_length)**i)
            
        # Populate arrays of neighbors accessible by point mutation, the intrinsic
        # probability of that error, and the number of steps away this neighbor is from
        # its parent
        num_seq = len(self.seq_list)
        self.neighbor_indexes = np.zeros((num_seq,num_neighbors),dtype=int)
        self.neighbor_fwd_prob_index[i,j,k] = np.ones((num_seq,num_neighbors,self.num_steps),dtype=int)
        self.neighbor_rev_prob_index[i,j,k] = np.ones((num_seq,num_neighbors,self.num_steps),dtype=int)

        # Populate np arrays that point from sequence index to neighbor indexes.  
        # Record indexes to the appropriate transition probability indexes as well.
        new_seq_list = self.seq_list[:]
        for i, seq in enumerate(seq_list):
            for j, neighbor in enumerate(possibleNeighbors(seq,self.num_steps,self.tm.alphabet)):

                neighbor_seq = neighbor[0]
                mutations = neighbor[1]

                # Record index of this neighbor in seq_list so we can point to it
                # in likelihood function 
                try:
                    self.neighbor_indexes[i,j] = self.seq_index_dict[neighbor_seq]


                # If the neighbor sequence was not actually seen, point to a 
                # dummy sequence.  There will be a dummy sequence for each 
                # class of possible mutations (e.g. a combination of A->T,P->Q
                # will have a different dummy sequence than A->T,R->S)
                except KeyError:

                    # Grab the set of all mutations
                    mut_key = []
                    for m in mutations:
                        mut_key.append((m[0],m[2]))
                    mut_key = tuple(mut_key)
                   
                    k = self.seq_index_dict[mut_key]
 
                    self.neighbor_indexes[i,j] = k
                    self.neighbor_indexes[k,xxx] = i

                    # Record it 
                    #try:

                    # but if this mutation class hasn't been seen yet, add it.
                    #except KeyError:

                    #    new_seq_list.append(mut_key)
                    #    self.seq_index_dict[mut_key] = len(new_seq_list) - 1
                    #    self.seq_counts.append(0)
                    #    self.neighbor_indexes[i,j] = self.seq_index_dict[mut_key]

                k = 0
                for m in mutations:
                    self.neighbor_fwd_prob_index[i,j,k] = self.tm.letter_to_index[m[0]]*self.tm.alphabet_length + letter_to_index[m[2]]
                    self.neighbor_rev_prob_index[i,j,k] = self.tm.letter_to_index[m[2]]*self.tm.alphabet_length + letter_to_index[m[1]]
                    k += 1

        self.seq_list = new_seq_list[:]
   
        self.transition_matrix = np.prod(err_prob 



        x = np.array(range(len(self.seq_list))
        y = np.array(self.seq_counts,dtype=float)

        y = param[x]
        y_out = param[x]*param[self.neighbor_fwd_prob_index[:,self.neighbor_indexes[x,:],:]]
        y_in  = param[x]*param[self.neighbor_rev_prob_index[:,self.neighbor_indexes[x,:],:]]
       
        y_obs = y - y_out + y_in 

        self.seq_counts*log(y_obs)
 
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
        """
    
        f = param[:self.size]
        Q = param[self.size]*f*self.transition_matrix
        return np.sum(np.log(f - np.sum(Q,1) + np.sum(Q,0))*self.seq_counts)


"""
class SequencingErrorProcessor(BaseProcessor):

    def process(self):

        pass

        # DO FIT FOR EACH ROUND INDEPENDENTLY...
        seq_list = count_dict.keys()
        seq_as_string = [converter(s) for s in seq_list]

        #round_list = [SOMEHOW POPULATE HERE]

        seq_dict = {}
        for i in range(len(seq_list)):
            seq_dict[seq_as_string[i]] = tuple([c for c in count_dict[seq_list[i]] if c != None])

    @property
    def data(self):

        pass

"""
