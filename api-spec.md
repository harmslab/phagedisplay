#phage display data format/api spec 
VERY ROUGH DRAFT

##Overview
This data format/api is designed for handling phage display experimental output
data in a consistent and efficient way.  It links raw fastq files all the way
through a processing pipeline.  It can also take simulated input (using its
pickled dictionary format) and treat it with the same pipeline we would use to
analyze an experiment.

##Some design details
* This pipeline is built on a generic `ProcessorBase` class.  This class holds 
actual experimental information and metadata in a consistent format.  Class 
instances can be nested.  
 * Each instance has a `data` property that spits out the interesting 
data from the class in a useful, programatically accessible way.  
 * Each instance has a `process` method which takes `**kwargs` and, potentially
the `data` property from the previous instance and processes/transforms the
data appropriatesly. 
 * Under the hood, everything is stored as a list/tuple in which the index indicates
the selection round.  0 is the initial library, 1 is after one round, 2 after two 
rounds etc.  Rounds in which no data were collected should have `NoneType`.  
 * The key format is a dictionary that keys sequences to counts/frequencies 
over rounds.  Each index in the tuple corresponds to a round (as above), with None
for missing entries.  The following dictionary shows how two peptides changed over 
rounds.  The first went from 25 in the initial library to 300 at round 2 and to
5000 and round 3. These dictionaries are written out as binary/python3 pickles.
```
    this_dict = {"ATGPRCTNKRDY":(25,None,300,5000),
                 "ACRNDPQENDWV":(400,None,30,0),
                 ...}
```

##Saving
Processor instances save out data in three ways.  
* `save.pickle` holds all information in the class.  It can only be
restored by creating a new MasterProcessor instance and loading this file via the
loadFile method.
* `info.json` holds all of the metadata in the class, dumped into a format that
can be read by people or other standard json readers.  Information is lost
relative to save.pickle, so this should not be used save and restore actual
Processor classes.
* Appropriate datasets. For example, FastqListProcessor stores gzipped fastq
files, while FastqToCountsProcessor stores pickled sequence/round dictionaries.  

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


##Example: load in fastq files and count how often each peptide occurs in each round

Create a master processor *m* that will store everything it does in the 
directory *test*.  
```
m = processors.MasterProcessor(expt_name="test")
m.create()
```

Create a new processor for loading in lists of fastq files.  This is then
appended to the master container, which will take charge of things like saving
the files, etc.  The final call of *m.process* takes a list of fastq_files which
the master processor will hand down to the appended subclass.  *process only 
operates on the most recently appended processor!*

```
p1 = processors.FastqListProcessor()
m.appendProcessor(p1)
m.process(list_of_fastq_to_load)
```

Create a new processor for converting those fastq files into something 
more useful (in this case, number of times each sequence is seen at each round).
The *m.process()* call doesn't take any arguments, as this class will draw from
the previous fastq file set.

```
p2 = processors.FastqToCountsProcessor()
m.appendProcessor(p2)
m.process()
```

