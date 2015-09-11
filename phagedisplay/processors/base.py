__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-05-16"

from . import blob
from .blob import Blob, DateBlob, FileBlob
import datetime, string, random, os, json, pickle

class BaseProcessor:
    """
    Base class for phage display processing steps.
    """

    def __init__(self,expt_name=None,date=None,description=None,ident=None):
        """
        """

        self._props = {}
        self._subprocessors = []
        self._do_not_write_to_json = []
        
        # Create a date stamp
        if not date:
            d = datetime.date.today()
            self.addProperty("date","{:d}-{:02d}-{:02d}".format(d.year,d.month,d.day),
                          DateBlob)
        else:
            self.addProperty("date",date,DateBlob)

        # Add a description
        self.addProperty("description",description)

        # Create a totally unique identifier string for this dataset
        if not ident:
            chars = string.ascii_letters + string.digits
            ident = "".join([random.choice(chars) for i in range(15)])
            self.addProperty("identifier",ident)
        else:
            self.addProperty("identifier",ident)
  
        # Add experiment name 
        if not expt_name:
            expt_name = ident
        self.addProperty("expt_name",expt_name) 


    def create(self,base_dir=None):
        """
        Create a directory in which the experiment will reside.
        """

        if base_dir:
            expt_name = os.path.join(base_dir,self.getProperty("expt_name"))
        else:
            expt_name = self.getProperty("expt_name")

        os.mkdir(expt_name)
        self.addProperty("expt_name",expt_name,FileBlob)       

        # Don't allow overwrite here, as this is when it is created (effectively 
        # a "save as" call)
        self.saveFile()

    def saveFile(self,filename=None,overwrite=False):
        """
        Write out the whole data structure tp a pickle.
        """

        # Write out pretty, human and other-programming-languages readable
        # json file.  
        self._writeJson(overwrite=overwrite)

        if not filename:
            filename = os.path.join(self.getProperty("expt_name"),"save.pickle")

        if not overwrite and os.path.exists(filename):
            raise IOError("{:s} file exists.".format(filename))

        f = open(filename,'wb')
        pickle.dump(self.__dict__,f)
        f.close() 

    def loadFile(self,filename):
        """
        Read a pickle file into an instance of the class.  (This will overwrite
        anything already stored in the instance that overlaps in __dict__, but
        keep the non-overlapping stuff).
        """

        f = open(filename,'rb')
        tmp_dict = pickle.load(f)
        f.close()

        self.__dict__.update(tmp_dict)

    def addProcessor(self,processor_class,**kwargs):
        """
        Add another processor within this one.  
        """

        if "expt_name" in kwargs:
            kwargs["expt_name"] = os.path.join(self.getProperty("expt_name"),
                                              kwargs["expt_name"])

        processor_class.create(self.getProperty("expt_name"))
        self._subprocessors.append(processor_class)
    
    def process(self,some_data):
        """
        Publicly accessible function to actually process the data.  This is 
        really only meaningful for subclasses that do more than store other
        processors.
        """

        pass 

            
    def getProperty(self,key):
        """
        Return the value associated with the blob in self._props[key].
        """
    
        try:
            return self._props[key].value
        except KeyError:
            return None
 
    def addProperty(self,key,value,custom_blob_class=None):
        """
        Add a property to the instance.  custom_blob_class can be used to define
        a custom property validator.
        """

        if custom_blob_class:
            self._props[key] = custom_blob_class(key,value)
        else:
            self._props[key] = Blob(key,value)


    def _writeJson(self,json_file="info.json",overwrite=False):
        """
        Create a json file and write out contents of self._props.
        """

        filename = os.path.join(self.getProperty("expt_name"),json_file)

        if not overwrite and os.path.exists(filename):
            raise IOError("file {:s} already exists".format(filename))

        # Create a dictionary of keys to write out
        write_dict = {}
        for k in self._props.keys():
            if k not in self._do_not_write_to_json:
                write_dict[k] = self._props[k].value
            else:
                write_dict[k] = "not written to json" 
      
        # Write out nested combination classes 
        out_subprocessors = []
        for s in self._subprocessors:
            s._writeJson(overwrite=True)
            out_subprocessors.append(s.getProperty("expt_name"))
        
        write_dict["_out_subprocessors"] = out_subprocessors       
 
        # Write out 
        f = open(filename,"w")
        json_out = json.dump(write_dict,f,sort_keys=True,indent=4) 
        f.close()       

    def _readJson(self,json_file="info.json"):
        """
        Read a json file into self._props.  Generally, this function should not
        be run, as self.saveFile/self.loadFile will preserve everything about
        this class instance (including custom class instances, etc.) while the
        json file will only have text references to the contents of those
        classes. In particular, every sub processor will have ProcessorParent
        class rather than the subclass it had on creation. Further, data in 
        self._do_not_write_to_json will be lost.  
        """
    
        f = open(json_file,"r")
        json_in = json.read(f)
        f.close()
        
        # Load in nested combination classes 
        out_subprocessors = json_in.pop("_out_subprocessors")
        for s in out_subprocessors:
            new_processor = ProcessorParent()
            new_processor.load(s)
            self._subprocessors.append(new_processor)
 
        for k in json_in:
            self.addProperty(k,json_in[k])

    @property
    def data(self):
        """
        Access this dataset in a programatically accessible way.
        """

        return None
        
    @property 
    def summary(self):
        """
        Return string summary of what's in this processor.
        """

        output = ["Instance of '{:s}' class with the properties:".format(self.__name__)]
        for p in self._props:
            output.append("    {:s} : {}".format(p,self._props[p].value))

        if len(self._subprocessors) > 0:
            output.append("\nClasses contained:")
            for s in self._subprocessors:
                output.append(s.summary())
            
        return "\n".join(output) 

