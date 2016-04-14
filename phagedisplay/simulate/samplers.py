__description__ = \
"""
Classes for doing sampling of pool instances.  
"""
__author__ = "Michael J. Harms, harmsm@gmail.com"
__date__ = "2015-04-25"

import numpy as np
import sys
from . import utility

class PoolIsEmptyError(Exception):
    """
    """

    pass

class BaseSampler(object):
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

        if sum(pool_instance.current_counts) == 0:
            raise PoolIsEmptyError
            
    
    def calcWeights(self,pool_instance,round_to_sample=-1):
        """
        Weight to place on a given sequence (for sampling with replacement).
        For most instances, this can just be frequency.  Subclasses may tweak
        this to include other processes.  
        """
       
        return (1.0*pool_instance.round_counts(round_to_sample))/(pool_instance.round_counts(round_to_sample).sum())
    
    def runExperiment(self,pool_instance,
                      sample_size,
                      round_to_sample=-1,
                      checkpoint=False,
                      append_to_current=True):
        """
        Run a sampling experiment.  
        
            Args: pool_instance (current pool)
                  sample_size (number of sequences to take forward)
                  round_to_sample (which round to sample. default of -1 takes last round)
                  checkpoint (whether we should note this as an important round)
                  append_to_current (if True, stick contents/counts on pool_instance. otherwise, return the values)

            Output: adds new round to pool_instance 
        
        """

        # Make sure it even makes sense to run an experiment
        self.poolSanityCheck(pool_instance)
        
        # If the sample_size is bigger than the number of sequences in the pool
        # and we don't allow replacement, just return the whole pool. 
        if sample_size > np.sum(pool_instance.round_counts(round_to_sample)) and not self.allow_replace:

            new_contents = pool_instance.round_contents(round_to_sample)
            new_counts = pool_instance.round_counts(round_to_sample)

        else:
            
            # If we allow replacement, do a simple weighted choice
            if self.allow_replace:

                x = self.calcWeights(pool_instance,round_to_sample)
                new_contents, new_counts = \
                    self._sample(pool_instance.round_contents(round_to_sample),
                                 weights=self.calcWeights(pool_instance,round_to_sample),
                                 sample_size=sample_size)
                
            # If we can't replace, create a possibilities array in which every
            # member of the pool is repeated as "counts" times.  Then calculate
            # weights accordingly.  
            else:
                possibilities = np.repeat(pool_instance.round_contents(round_to_sample),
                                          pool_instance.round_counts(round_to_sample))
                
                weights = self.calcWeights(pool_instance,round_to_sample)
                weights = np.repeat(weights,pool_instance.round_counts(round_to_sample))
                weights = weights/np.sum(weights)
                
                new_contents, new_counts = self._sample(possibilities,
                                                        weights=weights,
                                                        sample_size=sample_size)
      
        if append_to_current: 
            pool_instance.addNewStep(new_contents,new_counts,checkpoint)
        else:
            return new_contents, new_counts
        

class PCRAmplificationSampler(BaseSampler):
    """
    Class for simulating the amplification of phage that occurs via PCR.  
    """

    def __init__(self):
        
        self.allow_replace = True
    

class PhageAmplificationSampler(BaseSampler):
    """
    Class for simulating the amplification of phage that occurs after a binding
    round.
    """

    def __init__(self):
        
        self.allow_replace = True

class PipetteSampler(BaseSampler):
    """
    Class for simulating the downsampling that occurs when you take a small
    fraction of the total pool using a pipette.  
    """

    def __init__(self):
        
        self.allow_replace = False

    def calcWeights(self,pool_instance,round_to_sample=-1):
        """
        Since we're sampling without replacement, return even weights for every
        sequence.  Degeneracy will be taken care of in the sampling itself.  
        """
    
        return np.ones(pool_instance.round_counts(round_to_sample).size,dtype=float)/pool_instance.round_counts(round_to_sample).size
        
        
class BindingSampler(BaseSampler):
    """
    Class for simulating downsampling that occurs when the phage are placed on
    a plate.  
    """
    
    def __init__(self,conc_constant=1.0):
        """
        Args: 
            conc_constant scales affinities.
        """ 
 
        self.allow_replace = False
        self.conc_constant = conc_constant
    
    def calcWeights(self,pool_instance,round_to_sample=-1):
        """
        Weight by the affinity of the sequence.  Degeneracy will be taken care
        of since we sample without replacement.  
        """
    
        return pool_instance.round_affinities(round_to_sample)/np.sum(pool_instance.round_affinities(round_to_sample))
    
    def runExperiment(self,pool_instance,
                      sample_size,
                      round_to_sample=-1,
                      checkpoint=False,
                      append_to_current=True,
                      return_only_sample_size=False):
        """
        Run a sampling experiment.  
        
            Args: pool_instance (current pool)
                  sample_size (number of sequences to take forward)
                  round_to_sample (which round to sample. default of -1 takes last round)
                  checkpoint (whether we should note this as an important round)
                  append_to_current (if True, stick contents/counts on pool_instance. otherwise, return the values)
                  return_only_sample_size (if True, return sample size and do not sample)

            Output: adds new round to pool_instance 
        
        """

        # Make sure it even makes sense to run an experiment
        self.poolSanityCheck(pool_instance)
       
        # Grab the total number of sequences seen 
        total_counts = np.sum(pool_instance.round_counts(round_to_sample))

        weights = pool_instance.current_affinities*pool_instance.round_counts(round_to_sample)
        Q = np.sum(weights)

        # Calcualte probability that nothing binds at all (assumes that [M] is
        # much higher than Kx_average).
        p_nothing = self.conc_constant/(self.conc_constant + Q)

        # Scale sample size according to the probability of not binding
        sample_size = int(round(sample_size*(1-p_nothing[0]),0)) 

        if return_only_sample_size:
            return sample_size
  
        sys.stdout.write("# Number of proteins taken: {:d}\n".format(sample_size))
        sys.stdout.flush()
        # If the sample size remains bigger than the number of total sequences,
        # keep them all.  
        if sample_size >= sum(pool_instance.round_counts(round_to_sample)):
            new_contents = pool_instance.round_contents(round_to_sample)
            new_counts = pool_instance.round_counts(round_to_sample)

        # Otherwise, sample!
        else:
            weights = weights/Q

            # Normalize the binding probability to the number of times we will query
            # this probability in the final, expanded array
            weights = weights/pool_instance.round_counts(round_to_sample)

            # the number of times they occur.
            possibilities = np.repeat(pool_instance.round_contents(round_to_sample),
                                      pool_instance.round_counts(round_to_sample))
            weights = np.repeat(weights,pool_instance.round_counts(round_to_sample))

            # Sample                    
            new_contents, new_counts = self._sample(possibilities,
                                                    weights=weights,
                                                    sample_size=sample_size)

        if append_to_current:
            pool_instance.addNewStep(new_contents,new_counts,checkpoint)
        else:
            return new_contents, new_counts
    

class IlluminaSampler(BaseSampler):
    """
    """

    def __init__(self):
        
        self.allow_replace = True 
  

    def runExperiment(self,pool_instance,
                      sample_size,
                      round_to_sample=-1):
        """
        Sample from pool_instance at round_to_sample.  Return the contents
        and counts.  Do not append to pool.
        """
          
        ##### XXX ADD SEQUENCING ERRORS 
 
        return super(IlluminaSampler,self).runExperiment(pool_instance=pool_instance,
                                                         sample_size=sample_size,
                                                         round_to_sample=round_to_sample,
                                                         append_to_current=False)

