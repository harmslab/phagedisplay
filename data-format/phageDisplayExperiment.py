__description__ = \
"""
Data structure and functions for saving phage experimental and simulation
outputs in a common format.  For current format, see data-format-spec.txt in
this directory.
"""
__author__ = "Michael J. Harms"
__date__ = "2015-04-28"

DATA_FORMAT_VERSION = "0.1"
REQUIRED_KEYS = ["date","description","data_format_version","data_identifier"]
CORE_DIRECTORIES = ["fastqs","pickles"]

import os, shutil, json, random, string, datetime

class PhageDisplayExperiment:
    """
    Main data structure that holds phage display experiments, pointers to
    data files, and metadata. Function is to create a standardized directory
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
    
    def addFile(self,file_name,file_round,file_type="guess"):
        """
        Add a data file to the experiment.  
        
        Args:
            file_name: name of data file
            file_round: experimental round data file corresponds to (int >= 0)
            file_type: "pickle", "fastq" or "guess".  
        """

        # Make sure file_round is good
        if type(file_round) != int or (file_round < 0):
            err = "Round must be a positive integer."
            raise ValueError(err)

        # Grab file extension to guess
        if file_type == "guess":
            file_type = file_name.split(".")[-1]

        if file_type == "pickle":
            if not hasattr(self,"pickle_list"):
                self.pickle_list = []
            self.pickle_list.append((os.path.relpath(file_name),file_round))

        elif file_type == "fastq":
            if not hasattr(self,"fastq_list"):
                self.fastq_list = []
            self.fastq_list.append((os.path.relpath(file_name),file_round))

        else:
            err = "file type {:s} not recognized.".format(file_type)
            raise ValueError(err)

    def addProperty(self,key,value):
        """
        Add an arbitrary key/value property to the data structure.  

        Args: 
            key: property name (string)
            value: value (arbitrary data structure)
        """

        self._validateProperty(key,value)
        self._props[key] = value

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
        self._copyExternalFiles("fastq")
        self._copyExternalFiles("pickle")
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
            print("pickle_files: (filename, round)")
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

    def _copyExternalFiles(self,file_type):
        """
        Copy in any external files into the data directory.  In doing so, also
        update pickle_list and fastq_list so they have relative paths (base at
        self.dir_name).  

        Args:
            file_type: type of files being manipulated.  pickle or fastq
        """

        list_name = "{:s}_list".format(file_type)
        internal_dir_name = "{:s}s".format(file_type)

        if hasattr(self,list_name):
            for i, x in enumerate(self.__dict__[list_name]):

                # Try to copy in files.  If we die with SameFileError, we're 
                # copying to self and we can just ignore error.
                try:
                    new_name = os.path.join(internal_dir_name,
                                            os.path.split(x[0])[-1])
                    shutil.copy(x[0],os.path.join(self.dir_name,new_name))
                
                    # set relative path
                    self.__dict__[list_name][i] = (new_name,x[1])

                except shutil.SameFileError:
                    pass

    def _createInfoFile(self):
        """
        Create an info file and write out.
        """

        to_write = self._props.copy()

        if hasattr(self,"pickle_list"):
            for i, x in enumerate(self.pickle_list):
                to_write["pickle{:d}".format(i)] = x

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

            if k.startswith("pickle"):
                self._loadInternalFile("pickle",k,json_input[k])
            elif k.startswith("fastq"):
                self._loadInternalFile("fastq",k,json_input[k])
            else:
                self._props[k] = json_input[k]

    def _loadInternalFile(self,file_type,key,value):
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
