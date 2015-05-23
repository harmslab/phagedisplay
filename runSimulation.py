#!/usr/bin/env python
__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-05-21"
__usage__ = ""

import pickle
import simulate

e = simulate.StandardExperiment(sequence_length=7)
e.create(skew=3)
e.run()

pickle.dump(e.input_pool,open("input.pickle","wb"))
pickle.dump(e.actual_out,open("actual.pickle","wb"))
pickle.dump(e.illumina_out,open("illumina.pickle","wb"))
