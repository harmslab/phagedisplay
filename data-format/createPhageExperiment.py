#!/usr/bin/env python3
__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-28"

import argparse, sys, datetime

import processors

class XXXError(Exception):
    """
    """

    pass

def parseRoundsFile(rounds_file):
    """
    Grab the round numbers specified in a rounds_file, splitting on all spaces
    and line breaks.  Lines that start with "#" are ignored. 
    """

    f = open(rounds_file,'r')
    lines = f.readlines()
    f.close()

    # Skip commented lines
    lines = [l for l in lines if not l.startswith("#")]

    total_input = "  ".join(lines).split()
    round_numbers = [int(i) for i in total_input]
   
    return round_numbers   
 
   
def parseArbitraryKeyFile(arbitrary_key_file):
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
                raise XXXError(err)

            # Don't silently allow duplicate key/value pairs in the file
            if key in arb_key_dict:
                err = "key {:s} duplicated in key file {:s}\n".format(key,arbitrary_key_file)
                raise XXXError(err)

            arb_key_dict[key] = value

    return arb_key_dict

# Wrappers for class construction etc.
def createExperiment(expt_name,
                     description,
                     date,
                     fastq_files=[None],
                     pickle_file=None,
                     rounds_file=None,
                     arbitrary_key_file=None):
    """
    """

    # --------------------- fastq files ---------------------------------------
    # Create fastq_files input with None for skipped rounds
    if len([f for f in fastq_files if f != None]) > 0:

        if rounds_file:
            rounds = parseRoundsFile(rounds_file)
        else:
            rounds = range(len(fastq_files))

        new_fastq_files = [None for i in range(max(rounds)+1)]
        for i, f in enumerate(fastq_files):
            new_fastq_files[rounds[i]] = f
        
        fastq_files = new_fastq_files[:]

    # ---------------------- pickle files -------------------------------------
    if pickle_file:
        err = "pickle loading not yet implemented.  Sorry."
        raise XXXError(err)
   
    # ---------------------- other params -------------------------------------
    if date == "today":
        d = datetime.date.today()
        date = "{:04d}-{:02d}-{:02d}".format(d.year,d.month,d.day)

    arbitrary_keys = {}        
    if arbitrary_key_file:
        arbitrary_keys = parseArbitraryKeyFile(arbitrary_key_file)

    m = processors.MasterProcessor(expt_name=expt_name,
                                   date=date,
                                   description=description)

    # load in arbitrary meta data
    for k in arbitrary_keys:
        m.addProperty(k,arbitrary_keys[k])
   
    # Create the directory for the experiment 
    m.create()
    
    a = processors.FastqListProcessor(expt_name="raw-fastq-files")
    m.addProcessor(a)
    m.process(file_list=fastq_files)

    b = processors.FastqToCountsProcessor(expt_name="raw-counts")
    m.addProcessor(b)
    m.process()


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
    """

    if argv == None:
        argv = sys.argv[:]

    # Initialize argument parser
    parser = argparse.ArgumentParser(description='Create a phage experiment structure given input pickle and or fastq files.')

    # ----------------------- Positional arguments ----------------------------
    parser.add_argument('expt_name',
                        help="experiment name (will be directory name)")

    # ------------------------ Optional arguments  ----------------------------

    # You can specify pickle or fastq files, but not both
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
                           nargs=1,
                           help='processed/counted pickle file',
                           default=[None])

    # Sundry other arguments
    parser.add_argument('-d','--description',
                        dest="description",
                        nargs=1,
                        help="experiment description (string)",
                        default=["No description specified."])

    parser.add_argument('-t','--date',
                        dest="date",
                        nargs=1,
                        help="date in YYYY-MM-DD format (string)",
                        default=["today"])

    parser.add_argument('-r','--rounds',
                        dest="rounds_file",
                        nargs=1,
                        help="text file containing round numbers (if not specified, rounds are assumed sequential from 0).  Can only be specified for fastq files.",
                        action=FastqOnlyAction,
                        default=[None])

    parser.add_argument('-k','--keys',
                        dest="arbitrary_key_file",
                        nargs=1,
                        help="text file containing arbitrary key:value  pairs to put in data structure",
                        default=[None])

    # --------------------- Parse the command line ----------------------------
    args = parser.parse_args(argv[1:])

    p = createExperiment(expt_name=args.expt_name,
                         description=args.description[0],
                         date=args.date[0],
                         fastq_files=args.fastq_files,
                         pickle_file=args.pickle_file[0],
                         rounds_file=args.rounds_file[0],
                         arbitrary_key_file=args.arbitrary_key_file[0])


if __name__ == "__main__":
    main()
