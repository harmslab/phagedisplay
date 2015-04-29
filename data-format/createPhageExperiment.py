#!/usr/bin/env python3
__description__ = \
"""
Create a phage display experiment directory/file structure given a set of 
experimental inputs.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-28"

import argparse, sys, datetime
import phageDisplayExperiment

def createExperiment(out_dir,description,date,experiment_files):
    """
    """

    if date == "today":
        d = datetime.date.today()
        date = "{:04d}-{:02d}-{:02d}".format(d.year,d.month,d.day)

    p = phageDisplayExperiment.PhageDisplayExperiment(description=description,
                                                      date=date)

    fastq_counter = 0
    pickle_counter = 0
    for f in experiment_files:

        if f.split(".")[-1] == "pickle":
            p.addFile(f,pickle_counter)
            pickle_counter += 1
        
        if f.split(".")[-1] == "fastq":
            p.addFile(f,fastq_counter)
            fastq_counter += 1
    
    p.write(out_dir)
    p.print()

def main(argv=None):

    if argv == None:
        argv = sys.argv[:]

    parser = argparse.ArgumentParser(description='Create a phage experiment structure given input pickle and or fastq files.')

    # Positional arguments    
    parser.add_argument('out_dir',
                        help="output directory")
    parser.add_argument('input_files',
                        metavar='exp_data_file',
                        type=str,
                        nargs='+',
                        help='input fastq and/or pickle files, ordered by round')

    # Optional arguments
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

    # Parse the experiment
    args = parser.parse_args(argv[1:])
  
    out_dir = args.out_dir
    description = args.description[0]
    date = args.date[0]
    experiment_files = args.input_files

    createExperiment(out_dir=args.out_dir,
                     description=args.description[0],
                     date=args.date[0],
                     experiment_files=args.input_files)    

if __name__ == "__main__":
    main()
