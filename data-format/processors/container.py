__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-05-16"

import blob

class ContainerParent:
    """
    Base class for phage display processing steps.
    """

    def __init__(self):
        """
        """

        self._props = {}
        self._steps = []

    def create(date=None,description=None,ident=None):

        # Create a date stamp
        if not date:
            d = datetime.date.today()
            self._addProp("date","{:d}-{:02d}-{:02d}".format(d.year,d.month,d.day),blob.DateBlob)
        else:
            self._addProp("date",date,blob.DateBlob)

        # Add a description
        self._addProp("description",description)

        # Create a totally unique identifier string for this dataset
        if not ident:
            chars = string.ascii_letters + string.digits
            ident = "".join([random.choice(chars) for i in range(50)])
            self._addProp("identifier",ident)
        else:
            self._addProp("identifier",ident)

        # Create a base directory to store this stuff 
        if not base_dir:
            base_dir = ident
            os.mkdir(base_dir)

        self._addProp("base_dir",base_dir,blob.FileBlob)       

    def saveFile(self,filename=None):
        """
        """

        if not filename:
            filename = os.path.join(self.getProperty("base_dir"),"save.pickle")

        f = open(filename,'wb')
        pickle.dump(self.__dict__,f)
        f.close() 

    def loadFile(self,filename):
        """
        """

        f = open(filename,'rb')
        tmp_dict = pickle.load(f)
        f.close()

        self.__dict__.update(tmp_dict)

        #self._addProp(self.getProperty("base_dir"),base_dir,blob.FileBlob)
        #self._readJson(self,os.path.join(self.getProperty("base_dir"),"info.json"))
    
    def addStep(self,step_class,**kwargs):
        """
        """

        if "base_dir" in kwargs:
            kwargs["base_dir"] = os.path.join(self.getProperty("base_dir"),
                                              kwargs["base_dir"])

        new_step = step_class()
        new_step.create(**kwargs)
        self._steps.append(new_step)
    
    def process(self,some_data):
        """
        Publicly accessible function to actually process the data.
        """

        pass 

            
    def getProperty(self,key):
        """
        Return the value associated with the blob in self._props[key].
        """
    
        try:
            return self._props[k].value
        except KeyError:
            return None
 
    def _addProp(self,key,value,custom_blob_class=None):
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
            write_dict[k] = self._props[k].value
      
        # Write out nested combination classes 
        out_steps = []
        for s in self._steps:
            s._writeJson()
            out_steps.append(s.getProperty("base_dir"))
        
        write_dict["_out_steps"] = out_steps       
 
        # Write out 
        f = open(os.path.join(self.getProperty("base_dir"),json_file),"w")
        json_out = json.dump(to_write,f,sort_keys=True,indent=4) 
        f.close()       

    def _readJson(self,json_file="info.json"):
        """
        Read a json file into self._props
        """
    
        f = open(json_file,"r")
        json_in = json.read(f)
        f.close()
        
        # Load in nested combination classes 
        out_steps = json_in.pop("_out_steps")
        for s in out_steps:
            new_step = ContainerParent()
            new_step.load(s)
            self._steps.append(new_step)
 
        for k in json_in:
            self._addProp(k,json_in[k])

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

        if len(self._steps) > 0:
            output.append("\nClasses contained:")
            for s in self._steps:
                output.append(s.summary())
            
        return "\n".join(output) 

class TotalContainer(ContainerParent):
    """
    Holds onto an entire phage display experiment, with multiple container
    instances held therein.
    """
  
    @property
    def data(self):
        """
        Access this dataset in a programatically accessible way.
        """

        out = []
        for s in self._steps:
            out.append(s.data)

        return out
