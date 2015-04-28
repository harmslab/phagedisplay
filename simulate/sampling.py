__description__ = \
"""
Classes for doing sampling of pool instances.  
"""
__author__ = "Michael J. Harms, harmsm@gmail.com"
__date__ = "2015-04-25"

import numpy as np
import utility

class SamplerBaseClass(object):
    """
    Class for sampling pools.  
    
    Key attribute is "allow_replace."  If this is False, we are sampling from
    the whole pool.  If this is True, we are amplifying the pool.  
    """
    
    def __init__(self):
        """
        Initialize base class.  
        """
        
        self.allow_replace = True
        
    def _sample(self,possibilities,weights=np.array([]),sample_size=0):
        """
        Core sampling function.  
        """
        
        # If weights are specified, use them.
        if weights.size != 0:
            sampled = np.random.choice(possibilities,
                                       size=sample_size,
                                       replace=self.allow_replace,
                                       p=weights)
        else:
            sampled = np.random.choice(possibilities,
                                       size=sample_size,
                                       replace=self.allow_replace)
        
        # Collapse sample into a set of unique sequences with their counts
        return utility.uniqueCounter(sampled)
    
    def poolSanityCheck(self,pool_instance):
        """
        Make sure that it makes sense to sample this pool. 
        """
        
        # Has the pool been populated?
        if not pool_instance.pool_exists:
            err = "pool must be initialized prior to sampling.\n\n"
            raise PhageDisplaySimulatorError(err)
    
    
    def calcWeights(self,pool_instance):
        """
        Weight to place on a given sequence (for sampling with replacement).
        For most instances, this can just be frequency.  Subclasses may tweak
        this to include other processes.  
        """
        
        return pool_instance.current_counts/(pool_instance.current_counts.sum())
    
    def runExperiment(self,pool_instance,sample_size):
        """
        Run a sampling experiment.  
        
            Args: pool_instance (current pool)
                  sample_size (number of sequences to take forward)
        
            Returns: new_contents (new sequences after sampling)
                     new_counts (number of times each sequence is seen)
        
        """

        # Make sure it even makes sense to run an experiment
        self.poolSanityCheck(pool_instance)
        
        # If the sample_size is bigger than the number of sequences in the pool
        # and we don't allow replacement, just return the whole pool. 
        if sample_size > np.sum(pool_instance.current_counts) and not self.allow_replace:
            return pool_instance.current_contents, pool_instance.current_counts

        else:
            
            # If we allow replacement, do a simple weighted choice
            if self.allow_replace:
                new_contents, new_counts = \
                    self._sample(pool_instance.current_contents,
                                 weights=self.calcWeights(pool_instance),
                                 sample_size=sample_size)
                
            # If we can't replace, create a possibilities array in which every
            # member of the pool is repeated as "counts" times.  Then calculate
            # weights accordingly.  
            else:
                possibilities = np.repeat(pool_instance.current_contents,
                                          pool_instance.current_counts)
                
                weights = self.calcWeights(pool_instance)
                weights = np.repeat(weights,pool_instance.current_counts)
                weights = weights/np.sum(weights)
                
                new_contents, new_counts = self._sample(possibilities,
                                                        weights=weights,
                                                        sample_size=sample_size)
                
        return new_contents, new_counts
        

class PCRAmplificationSampler(SamplerBaseClass):
    """
    Class for simulating the amplification of phage that occurs via PCR.  
    """

    def __init__(self):
        
        self.allow_replace = True
    

class PhageAmplificationSampler(SamplerBaseClass):
    """
    Class for simulating the amplification of phage that occurs after a binding
    round.
    """

    def __init__(self):
        
        self.allow_replace = True

class PipetteSampler(SamplerBaseClass):
    """
    Class for simulating the downsampling that occurs when you take a small
    fraction of the total pool using a pipette.  
    """

    def __init__(self):
        
        self.allow_replace = False

    def calcWeights(self,pool_instance):
        """
        Since we're sampling without replacement, return even weights for every
        sequence.  Degeneracy will be taken care of in the sampling itself.  
        """
    
        return np.ones(pool_instance.current_counts.size,dtype=float)/pool_instance.current_counts.size
        
        
class BindingSampler(SamplerBaseClass):
    """
    Class for simulating downsampling that occurs when the phage are placed on
    a plate.  
    """
    
    def __init__(self):
        
        self.allow_replace = False
    
    def calcWeights(self,pool_instance):
        """
        Weight by the affinity of the sequence.  Degeneracy will be taken care
        of since we sample without replacement.  
        """
    
        return pool_instance.current_affinities/np.sum(pool_instance.current_affinities)
    
    def runExperiment(self,pool_instance,sample_size):
        """
        Run a sampling experiment.  
        
            Args: pool_instance (current pool)
                  sample_size (number of sequences to take forward)
        
            Returns: new_contents (new sequences after sampling)
                     new_counts (number of times each sequence is seen)
        
        """

        # Make sure it even makes sense to run an experiment
        self.poolSanityCheck(pool_instance)
       
        # Grab the total number of sequences seen 
        total_counts = np.sum(pool_instance.current_counts)
    
        # Calculate probability that each sequence binds to a single site.
        weights = pool_instance.current_affinities*pool_instance.current_counts
        Q = 1 + np.sum(weights)
        weights = weights/Q

        # Calcualte probability that nothing binds at all.
        p_nothing = 1/Q
      
        # Normalize the binding probability to the number of times we will query
        # this probability in the final, expanded array
        weights = weights/pool_instance.current_counts

        # Add a "nothing-bound" probability to the weights.  This will be added 
        # "total_counts" times, hence the normalization to total_counts.
        weights = np.append(p_nothing/total_counts,weights)
     
        # Add a "nothing bound" to the possible contents 
        contents = pool_instance.current_contents.copy()
        contents = np.append(-1,contents)

        # Add a total_counts "nothing bound" occurences to the counts 
        counts = pool_instance.current_counts.copy()
        counts = np.append(total_counts,counts)

        # Create possibilities and weights arrays arrays by repeating these instances
        # the number of times they occur.
        possibilities = np.repeat(contents,counts)
        weights = np.repeat(weights,counts)

        # Sample                    
        new_contents, new_counts = self._sample(possibilities,
                                                weights=weights,
                                                sample_size=sample_size)

        # toss empties if seen
        if new_contents[0] == -1:
            return new_contents[1:], new_counts[1:] 

        return new_contents, new_counts
    

class IlluminaRunSampler(SamplerBaseClass):
    pass
  

