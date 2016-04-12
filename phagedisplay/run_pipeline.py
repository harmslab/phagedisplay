#!/usr/bin/env python3
__description__ = \
"""
Run the processing pipeline on a set of input files.  (This may eventually be
the "master" script that we run on new input data).
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-28"

import argparse, sys, datetime, pickle, os

import phagedisplay
from phagedisplay import processors, simulate

class PhagePipelineError(Exception):
    """
    """

    pass


class Pipeline:
    """
    Class for running a series of common processes on phage display experimental
    and simulated input.
    """

    def __init__(self,
                 expt_name,
                 description="No description specified",
                 date="today",
                 arbitrary_key_file=None):
        """
        Initialize the class with common meta data.
        """

        if date == "today":
            d = datetime.date.today()
            date = "{:04d}-{:02d}-{:02d}".format(d.year,d.month,d.day)

        self.m = processors.MasterProcessor(expt_name=expt_name,
                                            date=date,
                                            description=description)

        self._parseArbitraryKeyFile(arbitrary_key_file)

        # Create the directory for the experiment 
        self.m.create()

    def loadInput(self):
        """
        Dummy in parent class.  Will be full-fledged method in childern.  After 
        running loadInput, the master processor -- self.m -- should have a 
        pickled dictionary of counts vs. round.
        """

        pass

    def finalize(self):
        """
        Given that the master processor has a dictionary of counts per round in 
        the last spot, do the rest of the processing business on it.
        """

        pass

    def _parseArbitraryKeyFile(self,
                               arbitrary_key_file=None):
        """
        Grab arbitrary keys/value pairs to be loaded into the main data set as
        metadata. 
        """
 
        arb_key_dict = {}
        if arbitrary_key_file:
            f = open(arbitrary_key_file,'r')
            lines = f.readlines()
            f.close()

            for l in lines:
                
                # Skip blank and commented lines 
                if l.startswith("#") or l.strip() == "":
                    continue
                    
                # Parse the line 
                try:
                    key = l.split(":")[0].strip()
                    value = ":".join(l.split(":")[1:]).strip()
                    if value == "":
                        raise IndexError

                except IndexError:
                    err = "Line:\n\n{:s}\n in {:s} does not have key:value format\n".format(l,arbitrary_key_file)
                    raise PhagePipelineError(err)

                # Don't silently allow duplicate key/value pairs in the file
                if key in arb_key_dict:
                    err = "key {:s} duplicated in key file {:s}\n".format(key,arbitrary_key_file)
                    raise PhagePipelineError(err)

                arb_key_dict[key] = value

            # load in arbitrary meta data
            for k in arb_key_dict:
                self.m.addProperty(k,arb_key_dict[k])

class FastqPipeline(Pipeline):
    """
    Load in fastq files and spit out a counts-style dictionary.
    """
    
    def loadInput(self,
                  fastq_files=[None],
                  rounds_file=None):
        """
        Load in fastq files (and maybe rounds files), process the raw sequencing
        input, and spit out some pretty sequence/counts-per-round dictionaries.
        """

        self.fastq_files = fastq_files
        self.rounds_file = rounds_file

        self._parseRoundsFile()

        # Create fastq_files input with None for skipped rounds

        new_fastq_files = [None for i in range(max(self.rounds)+1)]
        for i, f in enumerate(fastq_files):
            new_fastq_files[self.rounds[i]] = f
        
        fastq_files = new_fastq_files[:]
    
        # Load in the fastq files 
        a = processors.FastqListProcessor(expt_name="raw-fastq-files")
        self.m.addProcessor(a)
        self.m.process(file_list=fastq_files)
    
        # Count good sequences in the fastq files
        b = processors.FastqToCountsProcessor(expt_name="raw-counts")
        self.m.addProcessor(b)
        self.m.process()

        c = processors.RegressEnrichmentProcessor(expt_name="regression")
        self.m.addProcessor(c)
        self.m.process()
        
        d = processors.ClusterProcessor(expt_name="clustering")
        self.m.addProcessor(d)
        self.m.process()

    def _parseRoundsFile(self):
        """
        Grab the round numbers specified in a rounds_file, splitting on all spaces
        and line breaks.  Lines that start with "#" are ignored. 
        """

        if not self.rounds_file:
            self.rounds = range(len(self.fastq_files))
        else:
            f = open(self.rounds_file,'r')
            lines = f.readlines()
            f.close()

            # Skip commented lines
            lines = [l for l in lines if not l.startswith("#")]

            total_input = "  ".join(lines).split()
            self.rounds = [int(i) for i in total_input]
   

class PicklePipeline(Pipeline):
    """
    Rather simple pipeline, which loads in a pre-calculated counts-per-round 
    dictionary, previously stored as a pickle file.
    """  
 
    def loadInput(self, 
                  pickle_file=None):
        """
        Read in a counts-per-round dict from a dictionary.
        """

        # load in the pickle file
        a = processors.PickleDictProcessor(expt_name="input-pickle-file")
        self.m.addProcessor(a)
        self.m.process(pickle_file=pickle_file)

class SimulationPipeline(Pipeline):
    """
    Pipeline for loading in a simulated phage display experiment.
    """   
 
    def loadInput(self,
                  simulation_file=None,
                  all_samples=False):
        """
        Load in the simulation from "simualtion_file."  all_samples determines 
        whether we load in every sample from each round or instead use samples
        with simulated illumina sequencing.
        """
    
        # Load up the simulation 
        sim = simulate.StandardExperiment()
        sim.load(simulation_file)

        # Write out a temporary pickle file
        pickle_file = "tmp-sim-out.pickle" 
        if all_samples:
            pickle.dump(sim.actual_out,open(pickle_file,"wb"))
        else:
            pickle.dump(sim.illumina_out,open(pickle_file,"wb"))
  
        # Load in the pickle file 
        a = processors.PickleDictProcessor(expt_name="simulation-output")
        self.m.addProcessor(a)
        self.m.process(pickle_file=pickle_file)

        # delete the temporary file
        os.remove(pickle_file)
    
        c = processors.RegressEnrichmentProcessor(expt_name="regression")
        self.m.addProcessor(c)
        self.m.process()


# -----------------------------------------------------------------------------
# Stuff below here is for parsing the command line
# -----------------------------------------------------------------------------
   
class FastqOnlyAction(argparse.Action):
    """
    Subclass that makes sure that options associated with fastq lists only
    work if that option is specified.  Ripped from stackoverflow.

    https://stackoverflow.com/questions/11455218/python-argparse-enable-input-parameter-when-another-one-has-been-specified
    """

    def __call__(self,parser,namespace,values,option_string=None):

        was_fastq_set = getattr(namespace,'fastq_files',[None])

        if was_fastq_set == [None]:
            parser.error( "Option only comptaible with --fastq.")
        else:
            setattr(namespace,self.dest,values)


def main(argv=None):
    """
    Parse argv and run the pipeline.
    """

    if argv == None:
        argv = sys.argv[:]

    # Initialize argument parser
    parser = argparse.ArgumentParser(description='Create a phage experiment structure given input pickle and or fastq files.')

    # ----------------------- Positional arguments ----------------------------
    parser.add_argument('expt_name',
                        help="experiment name (will be directory name)")

    # ------------------------ Optional arguments  ----------------------------

    # You can specify pickle, fastq files or simulation file, but not together
    group_one = parser.add_mutually_exclusive_group()
    group_one.add_argument('-f','--fastq',
                           dest='fastq_files',
                           type=str,
                           nargs='*',
                           help='input fastq files, ordered by round',
                           default=[None])

    group_one.add_argument('-p','--pickle',
                           dest='pickle_file',
                           type=str,
                           help='processed/counted pickle file',
                           default=None)

    group_one.add_argument('-s','--simulation',
                           dest='simulation_file',
                           help='saved simulation output',
                           default=None)  

    # ------------------------ Globally useful optional arguments  ----------------------------

    parser.add_argument('-d','--description',
                        dest="description",
                        help="experiment description (string)",
                        default="No description specified.")

    parser.add_argument('-t','--date',
                        dest="date",
                        help="date in YYYY-MM-DD format (string)",
                        default="today")


    parser.add_argument('-k','--keys',
                        dest="arbitrary_key_file",
                        help="text file containing arbitrary key:value  pairs to put in data structure",
                        default=None)

    # ------------------------ fastq file specific optional args  ----------------------------

    parser.add_argument('-r','--rounds',
                        dest="rounds_file",
                        help="text file containing round numbers (if not specified, rounds are assumed sequential from 0).  Can only be specified for fastq files.",
                        action=FastqOnlyAction,
                        default=None)
    
    # ------------------------ simulation specific optional args  ----------------------------

    parser.add_argument('-a','--all-samples',
                        dest="all_samples",
                        help="use all samples from simulation, not just illumina-simulated samples.  (default: False) Only sensicle for a simulation input.",
                        action="store_true",
                        default=False)
    

    # --------------------- Parse the command line ----------------------------
    args = parser.parse_args(argv[1:])


    # List of fastq files specified
    if len([f for f in args.fastq_files if f]) > 0:
        
        p = FastqPipeline(expt_name=args.expt_name,
                          description=args.description,
                          date=args.date,
                          arbitrary_key_file=args.arbitrary_key_file)
        p.loadInput(fastq_files=args.fastq_files,
                    rounds_file=args.rounds_file)

    # pickle file specified 
    elif args.pickle_file:
    
        p = PicklePipeline(expt_name=args.expt_name,
                           description=args.description,
                           date=args.date,
                           arbitrary_key_file=args.arbitrary_key_file)
        p.loadInput(pickle_file=args.pickle_file)

    # simulation file specified
    elif args.simulation_file:

        p = SimulationPipeline(expt_name=args.expt_name,
                               description=args.description,
                               date=args.date,
                               arbitrary_key_file=args.arbitrary_key_file)

        p.loadInput(simulation_file=args.simulation_file,
                    all_samples=args.all_samples)

    # No file specified  
    else:
        
        p = Pipeline(expt_name=args.expt_name,
                     description=args.description,
                     date=args.date,
                     arbitrary_key_file=args.arbitrary_key_file)
        
        p.loadInput()

    # Do final analysis pipeline
    p.finalize()

if __name__ == "__main__":
    main()
