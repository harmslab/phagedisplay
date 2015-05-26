#!/usr/bin/env python3
__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-05-21"
__usage__ = ""

import argparse, pickle, sys, os
import simulate

class RunSimulationError(Exception):
    """
    General error class for this module.
    """

    pass

def checkFileExists(filename):
    """
    See file exists.  If so, raise error.
    """

    if os.path.exists(filename):
        err = "output file {:s} exists".format(filename)
        raise RunSimulationError(err)
    

def runSimulation(out_file="save.simulation",
                  num_rounds=3,
                  param_set="lwheeler00",
                  sequence_length=7,
                  k_distrib_skew=3,
                  affinity_max=1e6,
                  conc_const=-1):
    """
    Run a simulation, writing everything out to out_file.  This binary file
    can then be read back in, preserving every aspect of the simulation for
    later analysis.
    """

    # Check for the output file before running the simulation so we don't finish
    # evertyhing and realize we don't have anywhere to put it! 
    checkFileExists(out_file)

    e = simulate.StandardExperiment(sequence_length=sequence_length)
    e.create(param_set=simulate.parameters.__dict__[param_set],
             skew=k_distrib_skew,   
             affinity_max=affinity_max,
             specified_conc_const=conc_const)
    e.run(num_rounds=num_rounds)
    e.save(out_file)

def main(argv=None):
    """
    Parse command line and run the simulation.
    """

    if argv == None:
        argv = sys.argv[:]

    # Initialize argument parser
    parser = argparse.ArgumentParser(description='Create a phage experiment structure given input pickle and or fastq files.')

    # ----------------------- Positional arguments ----------------------------
    parser.add_argument('out_file',
                        help="output file for the simulation")

    # ------------------------ Optional arguments  ----------------------------
    parser.add_argument('-n','--num-rounds',
                        dest="num_rounds",
                        help="number of rounds to simulate",
                        type=int,
                        default=3)

    possible_keys = [s for s in simulate.parameters.__dict__ if not s.startswith("_")]
    help_message = "name of parameters to use.  Possibilities: \n{:s}".format("   \n".join(possible_keys))
    parser.add_argument('-p','--param',
                        dest="param_set",
                        help=help_message,  
                        type=str,
                        default="lwheeler00")

    parser.add_argument('-c','--conc-const',
                        dest="conc_const",
                        type=float,
                        help="concentration constant is a constant that scales affinities to place them on an absolute scale.  A negative value indicates that this should be optimized prior to running.",
                        default=-1)

    parser.add_argument('-l','--length',
                        dest="sequence_length",
                        help="length of the simulated peptide sequence",
                        type=int,
                        default=7)

    parser.add_argument('-s','--skew',
                        dest="k_distrib_skew",
                        help="amount of skew to add too the affinity distribution.  Higher numbers skew towards more weaker binders in the dataset.",
                        type=int,
                        default=3)

    parser.add_argument('-m','--max-k',
                        dest="k_max",
                        help="affinity of the maximum sequence",
                        type=float,
                        default=1e6)
    

    # --------------------- Parse the command line ----------------------------
    args = parser.parse_args(argv[1:])

    runSimulation(out_file=args.out_file,
                  num_rounds=args.num_rounds,
                  param_set=args.param_set,
                  sequence_length=args.sequence_length,
                  k_distrib_skew=args.k_distrib_skew,
                  affinity_max=args.k_max,
                  conc_const=args.conc_const)



if __name__ == "__main__":
    main()
