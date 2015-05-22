#!/usr/bin/env python
__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-05-21"
__usage__ = ""

import pickle, matplotlib
import matplotlib.pyplot as plt

from simulate import pool, samplers

def runMe(sequence_length,num_rounds=3,affinity_max=1e6,alphabet=pool.AMINO_ACIDS,
          library_size_scalar=2.8e-5,screen_scalar=1.5e-3,amplify_scalar=6.8e-3,
          conc_constant=1.5e18):
    """
    Simulate a basic binding experiment.  
    """
    
    scalar = len(alphabet)**(sequence_length)
    library_size = scalar*library_size_scalar
    
    pipet = samplers.PipetteSampler()
    screen = samplers.BindingSampler(conc_constant=conc_constant)
    amplify = samplers.PhageAmplificationSampler()

    p = pool.Pool(sequence_length=sequence_length,alphabet=alphabet)
    p.createScaledPool(library_size,affinity_max)

    try:
        print(0,p.current_counts.size,sum(p.current_counts))
        for i in range(num_rounds):
            pipet.runExperiment(p,library_size)
            screen.runExperiment(p,screen_scalar*scalar)
            amplify.runExperiment(p,amplify_scalar*scalar,checkpoint=True)

            print(i+1,p.current_counts.size,sum(p.current_counts))
    
    except pool.PoolIsEmptyError:
        pass

    return p
    
p = runMe(3,num_rounds=3,affinity_max=1e6)

f = open(out_file,'wb')
pickle.dump(p.round_count_dict,f)
f.close()


