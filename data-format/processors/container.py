__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-05-16"

import blob
import datetime, string, random, os, json, pickle

class ContainerParent:
    """
    Base class for phage display processing steps.
    """

    def __init__(self):
        """
        """

        self._props = {}
        self._subcontainers = []
        self._do_not_write_to_json = []

    def create(self,expt_name=None,date=None,description=None,ident=None):
        """
        Create a new experiment with some common meta data.  This will also
        create the directory in which the experiment will reside.
        """

        # Create a date stamp
        if not date:
            d = datetime.date.today()
            self._addProperty("date","{:d}-{:02d}-{:02d}".format(d.year,d.month,d.day),
                          blob.DateBlob)
        else:
            self._addProperty("date",date,blob.DateBlob)

        # Add a description
        self._addProperty("description",description)

        # Create a totally unique identifier string for this dataset
        if not ident:
            chars = string.ascii_letters + string.digits
            ident = "".join([random.choice(chars) for i in range(15)])
            self._addProperty("identifier",ident)
        else:
            self._addProperty("identifier",ident)

        # Create a base directory to store this stuff 
        if not expt_name:
            expt_name = ident

        os.mkdir(expt_name)
        self._addProperty("expt_name",expt_name,blob.FileBlob)       

    def saveFile(self,filename=None):
        """
        Write out the whole data structure tp a pickle.
        """

        # Write out pretty, human and other-programming-languages readable
        # json file.  
        self._writeJson()

        if not filename:
            filename = os.path.join(self.getProperty("expt_name"),"save.pickle")

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

    def addSubContainer(self,container_class,**kwargs):
        """
        Add another container within this one.  **kwargs will all be passed to 
        the "create" method of container_class.
        """

        if "expt_name" in kwargs:
            kwargs["expt_name"] = os.path.join(self.getProperty("expt_name"),
                                              kwargs["expt_name"])

        new_container = container_class()
        new_container.create(**kwargs)
        self._subcontainers.append(new_container)
    
    def processData(self,some_data):
        """
        Publicly accessible function to actually process the data.  This is 
        really only meaningful for subclasses that do more than store other
        containers.
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
 
    def _addProperty(self,key,value,custom_blob_class=None):
        """
        Add a property to the instance.  custom_blob_class can be used to define
        a custom property validator.
        """

        if custom_blob_class:
            self._props[key] = custom_blob_class(key,value)
        else:
            self._props[key] = blob.Blob(key,value)


    def _writeJson(self,json_file="info.json"):
        """
        Create a json file and write out contents of self._props.
        """

        # Create a dictionary of keys to write out
        write_dict = {}
        for k in self._props.keys():
            if k not in self._do_not_write_to_json:
                write_dict[k] = self._props[k].value
      
        # Write out nested combination classes 
        out_subcontainers = []
        for s in self._subcontainers:
            s._writeJson()
            out_subcontainers.append(s.getProperty("expt_name"))
        
        write_dict["_out_subcontainers"] = out_subcontainers       
 
        # Write out 
        f = open(os.path.join(self.getProperty("expt_name"),json_file),"w")
        json_out = json.dump(write_dict,f,sort_keys=True,indent=4) 
        f.close()       

    def _readJson(self,json_file="info.json"):
        """
        Read a json file into self._props.  Generally, this function should not
        be run, as self.saveFile/self.loadFile will preserve everything about
        this class instance (including custom class instances, etc.) while the
        json file will only have text references to the contents of those
        classes. (In particular, every sub container will have ContainerParent
        class rather than the subclass it had on creation).   
        """
    
        f = open(json_file,"r")
        json_in = json.read(f)
        f.close()
        
        # Load in nested combination classes 
        out_subcontainers = json_in.pop("_out_subcontainers")
        for s in out_subcontainers:
            new_container = ContainerParent()
            new_container.load(s)
            self._subcontainers.append(new_container)
 
        for k in json_in:
            self._addProperty(k,json_in[k])

    @property
    def data(self):
        """
        Access this dataset in a programatically accessible way.
        """

        return None
        
    @property 
    def summary(self):
        """
        Return string summary of what's in this container.
        """

        output = ["Instance of '{:s}' class with the properties:".format(self.__name__)]
        for p in self._props:
            output.append("    {:s} : {}".format(p,self._props[p].value))

        if len(self._subcontainers) > 0:
            output.append("\nClasses contained:")
            for s in self._subcontainers:
                output.append(s.summary())
            
        return "\n".join(output) 

