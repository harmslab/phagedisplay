__description__ = \
"""
Data structure and functions for saving phage experimental and simulation
outputs in a common format.  For current format, see data-format-spec.txt in
this directory.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-28"

DATA_FORMAT_VERSION = "0.11"
REQUIRED_KEYS = ["date","description","data_format_version","data_identifier"]
CORE_DIRECTORIES = ["fastqs"]  

import os, shutil, json, random, string, datetime, sys
import processFastq

class PhageDisplayExperimentError(Exception):
    """
    General error class for this module.
    """
    
    pass

class PhageDisplayExperiment:
    """
    Main data structure that holds phage display experiments, pointers to
    data files, and metadata. Purpose is to create a standardized directory
    and information file that can then be read in and out for analysis.
    """

    def __init__(self,**kwargs):
        """
        Initialize the data structure.  Any **kwargs are stored as arbitrary
        properties that will be written to the info.json file. 
        """
      
        self._required_keys = REQUIRED_KEYS[:]
        self._core_directories = CORE_DIRECTORIES[:]
        self._props = {}
        
        # Do automatic generation first, then read kwargs so the user settings
        # always wipe out the default
        d = datetime.date.today()
        self._props["date"] = "{:d}-{:02d}-{:02d}".format(d.year,d.month,d.day)
        self._props["description"] = "Replace this with a better description!"
        self._props["data_format_version"] = DATA_FORMAT_VERSION

        # Create a totally unique identifier string for this dataset
        chars = string.ascii_letters + string.digits
        ident = "".join([random.choice(chars) for i in range(50)])
        self._props["data_identifier"] = ident
        
        # Read in user-specified properties
        for k in kwargs.keys():
            self._validateProperty(k,kwargs[k])
            self._props[k] = kwargs[k]


    def addFastqFile(self,file_name,file_round):
        """
        Add a fastq file to the experiment.  
        
        Args:
            file_name: name of data file
            file_round: experimental round data file corresponds to (int >= 0)
        """

        # Make sure file_round is good
        if type(file_round) != int or (file_round < 0):
            err = "Round must be a positive integer."
            raise ValueError(err)

        if not hasattr(self,"fastq_list"):
            self.fastq_list = []
        self.fastq_list.append((os.path.relpath(file_name),file_round))

    def addPickleFile(self,pickle_key,file_name):
        """
        Add a pickle file to the experiment.
    
        Args:
            pickle_key: what to name this pickle file in the data structure
            file_name: name of data file
        """

        if hasattr(self,pickle_key):
            err = "Pickle file \"{:s}\" already specififed for attribute \"{:s}\".\n".format(self._props[pickle_key],pickle_key)
            raise PhageDisplayExperimentError(err)
        else:
            self._props[pickle_key] = os.path.relpath(file_name)

    def addProperty(self,key,value):
        """
        Add an arbitrary key/value property to the data structure.  

        Args: 
            key: property name (string)
            value: value (arbitrary data structure)
        """

        self._validateProperty(key,value)
        self._props[key] = value

    def fastqToPickle(self):
        """
        """

        good_seq = "good-seq.pickle"
        bad_seq = "bad-seq.pickle"

        filenames = [f[0] for f in self.fastq_list]
        round_numbers = [f[1] for f in self.fastq_list]

        p = processFastq.CountFastqSeq()
        p.processFastqFiles(filenames,round_numbers,good_seq,bad_seq)
    
        self.addPickleFile("good_pickle_file",good_seq)
        self.addPickleFile("bad_pickle_file",bad_seq)

    def write(self,directory):
        """
        Write out the whole data structure to "directory."
        """
        self.dir_name = directory

        # Don't wipe out any directory that actually has stuff in it.
        if os.path.exists(self.dir_name):
            if os.listdir(self.dir_name):
                err = "Directory ({:s}) already exists.\n".format(self.dir_name)
                raise ValueError(err)
        else:
            os.mkdir(self.dir_name)
          
        # Make core directories 
        for d in self._core_directories:
            os.mkdir(os.path.join(self.dir_name,d))
                 
        # Create info.json and copy in pickle/fastq files 
        self._moveExternalFiles("fastq")
        self._movePickleFile("good_pickle_file")
        self._movePickleFile("bad_pickle_file")
        self._createInfoFile()

    def read(self,directory):
        """
        Load a saved experiment into this data structure.

        Args:
            directory: directory containing a previously written experiment.
        """
        
        self.dir_name = directory
        self._readInfoFile()
        
        print("File read.")
        self.print()

    def print(self):

        if hasattr(self,"dir_name"):
            print("\nContents of",self.dir_name)
        else:
            print("\nContents")

        print("")
        for k in self._props.keys():
            print(k,":",self._props[k])

        if hasattr(self,"pickle_list"):
            print("")
            print("pickle_list: (filename, useful)")
            for p in self.pickle_list:
                print("   ",p)
        
        if hasattr(self,"fastq_list"):
            print("")
            print("fastq_files: (filename, round)")
            for f in self.fastq_list:
                print("   ",f)

        print("")

    def _validateProperty(self,key,value):
        """
        Validate a key/value pair.  It's a hack at the moment with validation
        of specific entries hard-coded.  These should really be something like
        handlers defined for each property, but ...time.
        """
        # Make sure the date looks like 2015-04-29 
        if key == "date":
            try: 
                datetime.datetime.strptime(value, "%Y-%m-%d")
            except (ValueError,TypeError):
                err = "Property for \"{}\" \"{}\" invalid.".format(key,value)
                raise ValueError(err)
      
        # Make sure the specified data format is older than or equal to the
        # current format.
        if key == "data_format_version": 
            if float(DATA_FORMAT_VERSION) < float(value):
                err = "Data format ({:s}) is newer than parser format ({:s})".format(value,DATA_FORMAT_VERSION)
                raise ValueError(err) 

    def _moveExternalFiles(self,file_type):
        """
        Copy in any external files into the data directory.  In doing so, also
        update pickle_list and fastq_list so they have relative paths (base at
        self.dir_name).  

        Args:
            file_type: type of files being manipulated. fastq only at the moment
        """

        list_name = "{:s}_list".format(file_type)
        internal_dir_name = "{:s}s".format(file_type)

        if hasattr(self,list_name):
            for i, x in enumerate(self.__dict__[list_name]):

                # Try to move in files.  If we die with SameFileError, we're 
                # moving to self and we can just ignore error.
                try:
                    new_name = os.path.join(internal_dir_name,
                                            os.path.split(x[0])[-1])
                    shutil.move(x[0],os.path.join(self.dir_name,new_name))
                
                    # set relative path
                    self.__dict__[list_name][i] = (new_name,x[1])

                except shutil.SameFileError:
                    pass

    def _movePickleFiles(self,pickle_key):
        """
        Copy in any external pickle files into the data directory.  In doing so,
        also good_pickle_file and bad_pickle_file so they have relative paths, 
        (base at self.dir_name)

        Args:
            pickle_key: key pointing to pickle file ("good_pickle_file" etc.)
        """

        internal_dir_name = ""

        if hasattr(self,pickle_key):
            pickle_file = self._props[pickle_key]

            # Try to move in files.  If we die with SameFileError, we're 
            # moving to self and we can just ignore error.
            try:
                new_name = os.path.join(internal_dir_name,
                                        os.path.split(pickle_file)[-1])
                shutil.move(pickle_file,os.path.join(self.dir_name,new_name))
            
                # set relative path
                self._props[pickle_key] = new_name

            except shutil.SameFileError:
                pass

    def _createInfoFile(self):
        """
        Create an info file and write out.
        """

        to_write = self._props.copy()

        if hasattr(self,"fastq_list"):
            for i, x in enumerate(self.fastq_list):
                to_write["fastq{:d}".format(i)] = x
     
        f = open(os.path.join(self.dir_name,"info.json"),"w")
        json_out = json.dump(to_write,f,sort_keys=True,indent=4) 
        f.close()       
 
    def _readInfoFile(self):
        """
        Read an info file in from an existing directory.
        Args:
            info_file: name of the json-formatted info file for this 
                       experiment. 
        """

        f = open(os.path.join(self.dir_name,"info.json"),'r')
        json_input = json.loads(f.read())
        f.close()
   
        for k in json_input.keys():

            if k.startswith("fastq"):
                self._checkInternalFile("fastq",k,json_input[k])
            else:
                self._props[k] = json_input[k]

    def _checkInternalFile(self,file_type,key,value):
        """
        Load an internal file from the experimental data structure, making
        sure the file exists.
        """
        
        list_name = "{:s}_list".format(file_type)
        dir_name = "{:s}s".format(file_type)

        if not hasattr(self,list_name):
            self.__dict__[list_name] = []

        filename = os.path.join(self.dir_name,value[0])
        if os.path.isfile(filename):
            self.__dict__[list_name].append((filename,value[1]))
        else:
            err = "\"{:s}\" does not exist.".format(filename)
            raise IOError(err)


# Wrappers for class construction etc.
def createExperiment(output,
                     description,
                     date,
                     fastq_files=[None],
                     pickle_file=None,
                     rounds_file=None,
                     arbitrary_key_file=None):
    """
    """

    if date == "today":
        d = datetime.date.today()
        date = "{:04d}-{:02d}-{:02d}".format(d.year,d.month,d.day)

    # Dump "None" entries
    fastq_files = [f for f in fastq_files if f]
    if pickle_file != None and len(fastq_files) != 0:
        err = "You cannot specify both a pickle file and fastq_files.\n"
        raise PhageDisplayExperimentError(err)

    # Read in the rounds_file, if specified
    round_numbers = range(len(fastq_files))
    if rounds_file != None:

        if pickle_file:
            sys.stderr.write("Pickle file specified.  Ignoring \"rounds\" file.\n")
            sys.stderr.flush()

        else: 
            f = open(rounds_file,'r')
            lines = f.readlines()
            f.close()
    
            # Skip commented lines
            lines = [l for l in lines if not l.startswith("#")]

            total_input = "  ".join(lines).split()
            round_numbers = [int(i) for i in total_input]
   
            num_exp_rounds = len(fastq_files) 
            if len(round_numbers) != num_exp_rounds:
                err = "Number of rounds must match number of experimental input files."
                raise PhageDisplayExperimentError(err)

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
                raise PhageDisplayExperimentError(err)

            # Don't silently allow duplicate key/value pairs in the file
            if key in arb_key_dict:
                err = "key {:s} duplicated in key file {:s}\n".format(key,arbitrary_key_file)
                raise PhageDisplayExperimentError(err)

            arb_key_dict[key] = value
    
    p = PhageDisplayExperiment(description=description,
                               date=date,
                               **arb_key_dict)

    for i, f in enumerate(fastq_files):
        p.addFastqFile(f,round_numbers[i])
    
    p.fastqToPickle()
    p.write(output)

    p.print()

def loadExperiment():
    """
    What should this even be?
    """

    return p.load

def writeExperiment():
    pass
