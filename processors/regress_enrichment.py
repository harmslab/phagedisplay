__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-01-20"

from . import BaseProcessor

import sys, time, copy
import pickle, scipy
import numpy as np

from scipy.optimize import minimize, curve_fit

class DegenerateSequences:
    """
    Class that holds onto a set of sequences that share the exact same count
    pattern.
    """

    def __init__(self,pattern,sequence):
        """
        Initialize with a pattern and first sequence. 
        """
        self._pattern = tuple(pattern)
        self._sequences = [sequence]

    def append(self,sequence):
        """
        Add a new sequence to the list.
        """
        self._sequences.append(sequence)
        
    @property
    def sequences(self):
        return self._sequences

    @property
    def pattern(self):
        return self._pattern 

  
class DegenerateSequenceCluster:
    """
    Class that holds a set of DegenerateSequences instances that were placed 
    into a cluster because their patterns were similar.  
    """

    def __init__(self):
        """
        Create a list to hold all of the DegenerateSequences instances.
        """

        self.degen_sequence_set = []
        self.patterns_loaded = False
        self.pattern_length = -1

    def addSequenceSet(self,s):
        """
        Add a DegenerateSequences instance to the list.
        """

        if not self.patterns_loaded:
            self.patterns_loaded = True
            self.pattern_length = len(s.pattern)
        else:
            if len(s.pattern) != self.pattern_length:
                err = "Patterns must all have same length."
                raise ValueError(err)

        self.degen_sequence_set.append(copy.copy(s))

    def calcMeanPattern(self):
        """
        """

        pattern_list = []
        for s in self.degen_sequence_set:
            pattern_list.append(s.pattern)

        unique = list(dict([(p,()) for p in pattern_list]).keys())

        self.mean_pattern = np.zeros((len(self.degen_sequence_set[0].pattern)),
                                     dtype=float)
        self.degeneracy = 0
        for d in self.degen_sequence_set:
            pattern = np.array(d.pattern)
            repeats = len(d.sequences)
            
            self.mean_pattern = self.mean_pattern + pattern*repeats
            self.degeneracy += repeats
            
        self.mean_pattern = self.mean_pattern/self.degeneracy

        return self.mean_pattern, self.degeneracy    


class CoarseGrainer:
    """
    Coarse-grain a set of raw counts into bins.  Bin size is based on identifying
    non-overlapping bins based on poisson counting error.  In other words:

         (error)  (error)
    bin1--------||--------bin2---------|

    Bin size increases as counts increase because, for a poisson counting process,
    the absolute error increases with inreased magnitude.  Error is estimated by
    sqrt(N). 
    """   

    def __init__(self,maximum_value,cg_width=1):
        """
        Initialize self.breaks to hold a set of breaks with which to classify 
        inputs.  Maximum value sets the highest value that can be coarse-
        grained.  Anything above this will be set to the maximum bin.  Scalar
        is multipled by sqrt(value) to set width of coarse graining.  
        """

        v = np.array(range(2,maximum_value+1),dtype=int)
        v_err_floor = v - cg_width*np.sqrt(v)

        self.breaks = [0,1]
        current_ceiling = self.breaks[-1] + cg_width*np.sqrt(self.breaks[-1])
        for i in range(len(v)):
            if current_ceiling < v_err_floor[i]:
                self.breaks.append(v[i])
                current_ceiling = v[i] + cg_width*np.sqrt(v[i])

    def coarseGrain(self,v):
        """
        Walk through the breaks in self.breaks and figure out where each index
        in a list lands in the breaks list.  Pretty hacked...
        """

        # Bins holds bin assinment for each position in v; v_indexes holds the
        # indexes in v that have *not* been already assigned.  0 is always 0 in
        # our scheme, so these start as assigned
        bins = [0 for j in range(len(v))]
        v_indexes = [x for x in range(len(v)) if v[x] != 0]
       
        # Go through the breaks... 
        for i in range(1,len(self.breaks)):

            # Go through the sites in v that are not yet binned
            for j in v_indexes:
                
                # If we're in the bin, call it
                if v[j] > self.breaks[i-1] and v[j] <= self.breaks[i]:
                    bins[j] = i
                    v_indexes.remove(j)

            if len(v_indexes) == 0:
                break

        # If we didn't call everything, they are higher than max; record them as
        # being in the highest bin
        if len(v_indexes) != 0:
            for j in v_indexes:
                bins[j] = len(self.breaks)-1

        return bins           

    def bins2counts(self,bins):
        """
        Convert a set of bins spit out by coarseGrain back into a set 
        of putative counts (coarse-grained, obviously). 
        """

        return [self.breaks[b] for b in bins] 


class FitModel:
    """
    Class for doing a non-linear regression on phage display output data.

    Here are a few ways to transform the data to condition the regression.

    A = ln(theta_x_0)
    B = ln(E_x)

    bounds on A: < 0
    bounds on B: None
    hard part is that sum(theta_0) <= 1.

    theta_x(i) = theta_x_0*E_x**i/sum(Q)
               = theta_x_0*exp(i*B)/sum(Q)
               = exp(A)*E_x**i/sum(Q)
               = exp(A + B*i)/sum(Q)

    ln(theta_x(i)) = A + B*i - ln(sum(Q))
                   = ln(theta_x_0) + B*i - ln(sum(Q))
                   = A + ln(E_x)*i - ln(sum(Q))
                   = ln(theta_x_0) + ln(E_x)*i - ln(sum(Q))
    """

    def __init__(self,patterns,degeneracy,rounds):
        """
        Using the data in patterns and degeneracy, create an observable set, 
        objective function, constraints, and bounds taht can then be minimized.
        """

        self.patterns = patterns
        self.degeneracy = degeneracy 
        self.rounds = np.array(rounds,dtype=int)
        self.num_rounds = len(self.rounds)
        
        self.num_patterns = len(patterns)
        self.log_degeneracy = np.log(self.degeneracy)

        # Create observable matrix (thetas x rounds)
        self.y_obs = patterns[:]

        # Normalize observations so the frequency at each round adds up to 1.
        for j in range(self.num_rounds):
            Q = self.y_obs[:,j]*self.degeneracy
            self.y_obs[:,j] = Q/np.sum(Q)

        # Initialize temporary arrays used in objective function
        self.y_calc = np.zeros((self.num_patterns,self.num_rounds),
                               dtype=float)
        self.p = np.zeros((self.num_patterns),dtype=float)

        # Populate initial guesses
        # 1) Try to fit a simple single-site model to the data.  This is the
        #    first guess.
        # 2) If that doesn't converge, assign the initial conc to the frequency
        #    at obs0 and  K to 1.0
        print("Generating initial parameter guesses...",end="")
        sys.stdout.flush()
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

        print("Done.")
        sys.stdout.flush()

        # Create bound list.  conc must be between 0 and 1, K must be positive
        self.bounds = [(None,None) for i in range(self.num_patterns)]   
        self.bounds.extend([(None,None) for i in range(self.num_patterns)])

        #self.bounds[0] = (-1.0,-1.0)
        #self.param_guess[0] = -1.0
        self.bounds[self.num_patterns] = (-1.0,-1.0)
        self.param_guess[self.num_patterns] = -1.0
 
        # Constrain the sum of the initial conc parameters to be 1.
        #self.constraints = ({'type': 'eq',
        #                     'fun': lambda x:  1 - sum(x[:int(len(x)/2)])})
        

    def _objective(self,param):
        """
        Objective function to minimize.  Goes like:

        For each pattern i, the observed frequency (theta) at round j is:
 
            theta[i,j] = (theta[i,0]*(beta[i]**j))/sum_over_i_at_j

        Objective function is the rmsd between calculated theta and observed theta
        over all rounds.

        param = [all_thetas...., all_Ks...]
        """

        # Initial state 
        self.p = param[:self.num_patterns]

        # For every round, apply the beta term
        for j, n in enumerate(self.rounds):
            if n != 0:
                self.p = self.p*(param[self.num_patterns:])
            self.y_calc[:,j] = self.p/sum(self.p)

        return np.sum(np.power(self.y_calc-self.y_obs,2))


    def _objective2(self,param):
        """
        Objective function to minimize.  Goes like:

        ln(theta_x(i)) = A + B*i - ln(sum(Q))
        exp(A + B*i)/sum(Q)

        Objective function is the rmsd between calculated theta and observed theta
        over all rounds.

        param = [all_thetas...., all_Ks...]
        """

        for j, n in enumerate(self.rounds):
            self.p = np.exp(self.log_degeneracy + param[:self.num_patterns] + param[self.num_patterns:]*n)
            self.y_calc[:,j] = self.p/np.sum(self.p)

        return np.sum(np.power(self.y_calc-self.y_obs,2))


    def runRegression(self,maxiter=100000000000):
        """
        Run a regression of the objective function given our bound, constraints,
        etc.
        """

        self.start_time = time.time()

        self.fit_result = minimize(fun=self._objective2,x0=self.param_guess,
                                   bounds=self.bounds,
                                   #constraints=self.constraints,
                                   options={"maxiter":maxiter,"maxfun":maxiter})
        self.end_time = time.time()

        print(time.asctime())
        print(self.fit_result.status,self.fit_result.fun,self.fit_result.message)
        print(time.asctime())

    def returnParam(self):
        """
        """

        initial_theta = self.fit_result.x[:self.num_patterns]
        K_values = self.fit_result.x[self.num_patterns:]

        out = np.zeros((len(K_values),3),dtype=float)
        for i in range(3):
            out[:,i] = np.exp(self.log_degeneracy + initial_theta + K_values*(i+1))
            out[:,i] = out[:,i]/sum(out[:,i])

        return initial_theta, K_values, out   


    def _indepFit(self,x,conc=1.0,K=1.0):
    
        return conc*(K**x)


class RegressEnrichmentProcessor(BaseProcessor):
    """
    """

    def process(self,count_dict,cg_width=1.0,minimum_times_seen=4):
        """
        """

        self.count_dict = count_dict
        self.cg_width = cg_width
        self.minimum_times_seen = minimum_times_seen

        rounds, patterns = self._loadData()
    
        self.rounds = rounds
        self.patterns = patterns

        fit_pattern = np.zeros((len(self.patterns),self.patterns[0].pattern_length),dtype=float)
        degeneracy = np.zeros((len(self.patterns)),dtype=int)
        for i, c in enumerate(self.patterns):
            if len(c.degen_sequence_set) != 0:
                mean_pattern, degen = c.calcMeanPattern()
                fit_pattern[i,:] = mean_pattern
                degeneracy[i] = degen

        m = FitModel(fit_pattern,degeneracy,rounds)
        m.runRegression()

        theta, K, calc_values = m.returnParam()

        totals = np.zeros((3),dtype=float)
        for i, c in enumerate(self.patterns):
            if len(c.degen_sequence_set) != 0:
                for j in range(len(c.degen_sequence_set)):
                    for k in range(len(c.degen_sequence_set[j].sequences)):
                        totals += c.degen_sequence_set[j].pattern

        for i, c in enumerate(self.patterns):
            if len(c.degen_sequence_set) != 0:
                for j in range(len(c.degen_sequence_set)):
                    for k in range(len(c.degen_sequence_set[j].sequences)):
                        print(c.degen_sequence_set[j].sequences[k],
                              theta[i]/degeneracy[i],
                              K[i],
                              c.degen_sequence_set[j].pattern/totals,
                              calc_values[i,:]/degeneracy[i])

    def _loadData(self):
        """
        """
        
        all_sequences = list(self.count_dict.keys())

        # Create a list, (degen_sequence_list) that has a list of all sequences that
        # share exactly the same count pattern.  Each unique pattern/sequence set is
        # stored in an instance of the DegenerateSequences class.
        self.degen_sequence_dict = {}
        rounds = [i for i, v in enumerate(self.count_dict[all_sequences[0]])
                  if v != None]
        num_rounds = len(rounds)
        for s in all_sequences:

            # Get rid of missing data
            pattern = tuple([v for v in self.count_dict[s] if v != None])
            if len(pattern) != num_rounds:
                err = "All sequences must have same rounds observed."
                raise ValueError(err)

            if sum(pattern) < self.minimum_times_seen:
                continue
            try:
                self.degen_sequence_dict[pattern].append(s)
            except KeyError:
                self.degen_sequence_dict[pattern] = DegenerateSequences(pattern,s)

        degen_sequence_list = []
        for k in self.degen_sequence_dict.keys():
            degen_sequence_list.append(copy.copy(self.degen_sequence_dict[k]))

        num_unique_patterns = len(degen_sequence_list)
        print("Found %i unique patterns for %i clones" % (num_unique_patterns,
                                                          len(all_sequences)))


        # Create a numpy array of all unique patterns
        unique_patterns = np.zeros((num_unique_patterns,num_rounds),dtype=int)
        for i in range(num_unique_patterns):
            unique_patterns[i,:] = np.array(degen_sequence_list[i].pattern)

        # Coarse grain these unique patterns into similar patterns
        cg = CoarseGrainer(np.max(unique_patterns),self.cg_width)
        similar_pattern_dict = {}
        for i in range(num_unique_patterns):

            key = tuple(cg.coarseGrain(unique_patterns[i,:]))
            try:
                similar_pattern_dict[key].append(degen_sequence_list[i])
            except KeyError:
                similar_pattern_dict[key] = [degen_sequence_list[i]]

        # Create a final set of clusters (made of SequenceCluster instance)
        final_clusters = [DegenerateSequenceCluster()
                          for i in range(len(similar_pattern_dict))]
        for i, key in enumerate(similar_pattern_dict.keys()):
            for cluster in similar_pattern_dict[key]:
                final_clusters[i].addSequenceSet(cluster)

        print("Collapsed to %i similar patterns." % (len(final_clusters)))

        return rounds, final_clusters

 

