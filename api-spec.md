Title: Phage display experiment data structure spec.
Version: 0.11
Date: 2015-04-29
Author: Michael J. Harms

##Overview
This pipeline is built on a generic "processor" class which stores data for
each processing step as well as meta-data in the form of "blobs." 

###Example
```
m = processors.MasterProcessor(expt_name="test")
m.create()

p1 = processors.FastqListProcessor()
m.appendProcessor(p1)
m.process(list_of_fastq_to_load)

p2 = processors.FastqToCountsProcessor()
m.appendProcessor(p2)
m.process()

...
```
 
###Breakdown
```
m = processors.MasterProcessor(expt_name="test")
m.create()
```

This creates a master processor *m* that will store everything it does in the 
directory *test*.  

```
p1 = processors.FastqListProcessor()
m.appendProcessor(p1)
m.process(list_of_fastq_to_load)
```

This creates a new processor for loading in lists of fastq files.  This is then
appended to the master container, which will take charge of things like saving
the files, etc.  The final call of *m.process* takes a list of fastq_files which
the master processor will hand down to the appended subclass.  *process only 
operates on the most recently appended processor!*

```
p2 = processors.FastqToCountsProcessor()
m.appendProcessor(p2)
m.process()
```

This creates a new processor for converting those fastq files into something 
more useful (in this case, number of times each sequence is seen at each round).
The *m.process()* call doesn't take any arguments, as this class will draw from
the previous fastq file set.


##Saving
Processor instances save out data in three ways.  
* "save.pickle" holds all information in the class.  It can only be
restored by creating a new MasterProcessor instance and loading this file via the
loadFile method.
* "info.json" holds all of the metadata in the class, dumped into a format that
can be read by people or other standard json readers.  Information is lost
relative to save.pickle, so this should not be used save and restore actual
Processor classes.
* Appropriate datasets. For example, FastqListProcessor stores gzipped fastq
files.  

## Stored file types:
### pickled sequence files
A pickled dictionary (python3, binary) containing sequences as keys and tuples
of counts for each run.  First entry in tuple is round 0.  Missing rounds have 
None.  

### fastq
fastq files come to us from the sequencing facility.  The "round" attribute
specifies where this sample comes from during the selection rounds.  These are
zipped via the gzip library.  

```
An entry in one of these files has the form:

@NS500451:25:H3K2MBGXX:1:11101:11677:1043 1:N:0:CTCTCTAC+GCGATCTA
GCGTGGCAGATTCCTTATAATGCGTATGNTNNTNNTGNNNGNNNTTNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
+
AAAAAFFFFFFFFFFFFFAF.FFFFFFF#F##F##FF###<###FF#############################
```

We care (mostly) about the second line, which has the actual sequence.

