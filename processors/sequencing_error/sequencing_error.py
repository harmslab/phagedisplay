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

    def __init__(self,count_dict,num_steps,transition_matrix,err_guess=0.01):
        """
        Do other stuff
        
        Create a matrix for transitions between each possible sequence neigbor
        with a probability given by the types of transitions between each
        neighbor.
        """

        self.num_steps = num_steps
        self.tm = transition_matrix
        self.err_guess = err_guess

        # Convert count_dictionary to a list, sorted from highest to lowest counts
        # Toss any sequences with 0 counts or None
        tmp_seq_list = list((v[1],v[0]) for v in count_dict.items())
        tmp_seq_list = [s for s in tmp_seq_list if s[0] < 1 or s[0] == None]
        tmp_seq_list.sort(reverse=True)

        # Create core mapping arrays/dictionaries to link sequences back to their
        # indexes, counts, etc.
        self.num_real_seq = len(tmp_seq_list)
        self.seq_list =   [s[1] for s in tmp_seq_list]
        self.seq_counts = [s[0] for s in tmp_seq_list]
        self.seq_index_dict = dict((s,i) for i,s in enumerate(self.seq_list))

        # Create lists that will store all neighbors and the probability of reaching
        # them given the transition matrix
        self.neighbor_lists = [[] for s in self.seq_list]
        self.neighbor_prob = [[] for s in self.seq_list]
        self.neighbor_num_steps = [[] for s in self.seq_list]

        # Go through every sequence and calculate possible transitions to every 
        # possible neighbor
        tmp_seq_list = self.seq_list[:]
        for i, seq in enumerate(tmp_seq_list):
            for neighbor in possibleNeighbors(seq,self.num_steps,self.tm.alphabet):

                # Sequence and mutations required to get there
                neighbor_seq = neighbor[0]
                mutations = neighbor[1]

                # Grab the index of the neighbor
                try:
                    j = self.seq_index_dict[neighbor_seq]

                # If that neighbor hasn't been seen yet, create it.
                except KeyError:

                    self.seq_list.append(neighbor_seq)
                    self.seq_counts.append(0)
                    self.seq_index_dict[neighbor_seq] = len(self.seq_list)-1
                    self.neighbor_lists.append([])
                    self.neighbor_prob.append([])
                    self.neighbor_num_steps.append([])
            
                    j = self.seq_index_dict[neighbor_seq]

                # Record neighbor adjacency both ways
                self.neighbor_lists[i].append(j)
                self.neighbor_lists[j].append(i)

                # Record probability of errors going each direction from the
                # transition matrix
                self.neighbor_prob[i].append(1.0)
                self.neighbor_prob[j].append(1.0)
                for m in mutations:
                    a = self.tm.letter_to_index[m[0]]
                    b = self.tm.letter_to_index[m[2]]

                    self.neighbor_prob[i][-1] *= self.tm.transition_matrix[a,b]
                    self.neighbor_prob[j][-1] *= self.tm.transition_matrix[b,a]

                self.neighbor_num_steps[i].append(len(mutations))
                self.neighbor_num_steps[j].append(len(mutations))

        # Construct a total transition matrix for each i-->j transition
        self.total_matrix = np.zeros((len(self.seq_list),len(self.seq_list)),dtype=float)
        self.num_steps_matrix = np.zeros((len(self.seq_list),len(self.seq_list)),dtype=int)
        for i in range(len(self.seq_list)):
            self.total_matrix[self.neighbor_lists[i],i] = self.neighbor_prob[i]
            self.num_steps_matrix[self.neighbor_lists[i],i] = self.neighbor_num_steps[i]
      
        # This information has now been matrix-ified.  Delete.   
        del self.neighbor_lists
        del self.neighbor_prob 

        self.param_guesses = np.zeros((len(self.seq_list)+1),dtype=float)
        self.param_guesses[range(len(self.seq_list))] = self.seq_counts
        self.param_guesses = self.param_guesses/np.sum(self.param_guesses)
        self.param_guesses[-1] = self.err_guess

    def runRegression(self,maxiter=100000000000):
        """
        Run a regression of the objective function.
        """

        self.fit_result = minimize(fun=self._objective,
                                   x0=self.param_guess,
                                   options={"maxiter":maxiter,"maxfun":maxiter})

    def returnParam(self):
        """
        Return the fit parameters. 
            out is a list of tuples that encode sequence, estimated freq, and
                how many times this sequence was actually seen.
            err_rate is the estimated global error rate applied to the 
                transition matrix.
        """

        frequencies = self.fit_result.x[:-1]
        err_rate = self.fit_result.x[-1]

        out = []
        for i in range(len(self.seq_list)):
            out.append((self.seq_list[i],frequencies[i],self.seq_counts[i]))
  
        return out, err_rate

    def _objective(self,param):
        """
        Calculate the -log likelihood of the experimentally observed sequence
        counts given our model of mutational neighbors.
        """
        
        # Frequencies 
        f = param[:self.num_total_seq]

        # Q is the frequency of each sequence times the estimated error rate 
        # raised to the number of steps made, times the total error transition
        # matrix.
        Q = (param[-1]**self.num_steps_matrix)*f*self.total_matrix

        # Return -log(likelihood)
        return -1*np.sum(np.log(f - np.sum(Q,1) + np.sum(Q,0))*self.seq_counts)


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
