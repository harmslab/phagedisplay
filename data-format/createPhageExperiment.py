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

def main(argv=None):

    if argv == None:
        argv = sys.argv[:]

    parser = argparse.ArgumentParser(description='Create a phage experiment structure given input pickle and or fastq files.')

    # Positional arguments    
    parser.add_argument('output',
                        help="output directory")

    # Optional arguments
    parser.add_argument('-f','--fastq',
                        dest='fastq_files',
                        type=str,
                        nargs='*',
                        help='input fastq files, ordered by round',
                        default=[None])
    
    parser.add_argument('-p','--pickle',
                        dest='pickle_file',
                        type=str,
                        nargs=1,
                        help='processed/counted pickle file',
                        default=[None])

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
                        help="text file containing round numbers (if not specified, rounds are assumed sequential from 0)",
                        default=[None])

    parser.add_argument('-k','--keys',
                        dest="arbitrary_key_file",
                        nargs=1,
                        help="text file containing arbitrary key:value  pairs to put in data structure",
                        default=[None])

    # Parse the experiment
    args = parser.parse_args(argv[1:])

    if args.pickle_file[0] != None:
        if args.rounds_file[0] != None:
            err = "You cannot supply both a rounds file and a pickle file.\n"
            #raise phageDisplayExperiment.PhageDisplayExperimentError(err)
        if args.fastq_files[0] != None:
            err = "You cannot supply both fastq files and a pickle file.\n"
            #raise phageDisplayExperiment.PhageDisplayExperimentError(err)
    
    p = phageDisplayExperiment.createExperiment(output=args.output,
                                                description=args.description[0],
                                                date=args.date[0],
                                                fastq_files=args.fastq_files,
                                                pickle_file=args.pickle_file[0],
                                                rounds_file=args.rounds_file[0],
                                                arbitrary_key_file=args.arbitrary_key_file[0])


if __name__ == "__main__":
    main()
