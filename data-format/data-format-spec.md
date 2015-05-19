Title: Phage display experiment data structure spec.
Version: 0.11
Date: 2015-04-29
Author: Michael J. Harms

## Rewriting
Current format is trying to be too many things to too many people.  Going to
shift to a new strategy:

"""
Container instance, which is made of:
    ProcessStep instances, which are made of:
        Blob instances for different data types

c = Container()
SomeProcessStep(c,process_args)
SomeOtherProcessStep(c,other_process_args)
"""       

The container will have a directory with an info.json file.  Each process step
will then have its own info.json file with metadata/parameters, as well as any
output files.  (The input files will be in the previous step).  For fastq and
pickle imports, the original file will be outside of the data structure, however
these operations copy in the file rather than modify it, so the "chain of 
custody" is preserved.

 

## Overview 
This file/data structure will be generated and populated by a script given input
fastq and/or pickle files.  This can be done from either experimental or
simulated data. 

#### Example:
"""
    experiment-directory/
        info.json
        fastqs/
            round0.fastq
            round1.fastq
        good-reads.pickle
        bad-reads.pickle
"""

### experiment-directory
This directory can have any name.

### info.json
info.json contains meta data about the experiment.  It requires entries date, 
description, data_format_version, and data_identifier.  It can have entries
specifying local fastq and pickle file locations (optional, but without them the
experiment is empty). It can also store arbitrary key names that the python data
structure will be smart enough to read.  Examples of arbitrary keys might be 
round barcodes, person who worked at the bench, exact experimental conditions,
etc.

###Example:
"""
info.json
    {
        "date":"2015-04-30",
        "description":"description here",
        "data_format_version":"0.1",
        "data_identifier":"long_unique_automatically_generated_random_string",

        "arbitrary_key":"has a value", 

        "good_pickle_file":"good-reads.pickle",
        "bad_pickle_file":"bad-reads.pickle",
 
        "fastq0":{
                  "file":"fastq/round0.fastq",
                  "round":0
                  },
        "fastq1":{
                  "file":"fastq/round5.fastq",
                  "round":5
                  }
    }
"""

### good_pickle_file
A pickled dictionary (python3, binary) containing sequences as keys and tuples
of counts for each run.  First entry in tuple is round 0.  Missing rounds have 
None.  

### fastq
fastq files come to us from the sequencing facility.  The "round" attribute
specifies where this sample comes from during the selection rounds.

####Example:
"""
An entry in one of these files has the form:

@NS500451:25:H3K2MBGXX:1:11101:11677:1043 1:N:0:CTCTCTAC+GCGATCTA
GCGTGGCAGATTCCTTATAATGCGTATGNTNNTNNTGNNNGNNNTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+
AAAAAFFFFFFFFFFFFFAF.FFFFFFF#F##F##FF###<###FF#############################

We care (mostly) about the second line, which has the actual sequence.
"""

