__description__ = \
"""
"""
__author__ = "Michael J. Harms"
__date__ = "2015-05-16"


class ProcessorBase:
    """
    Base class for phage display processing steps.
    """

    def __init__(self):
        """
        """

        self._props = {}


    def create(date=None,description=None,ident=None):

        if not date:
            d = datetime.date.today()
            self._addProp("date","{:d}-{:02d}-{:02d}".format(d.year,d.month,d.day),blob.DateBlob)
        else:
            self._addProp("date",date,blob.DateBlob)

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

    def load(base_dir):
        """
        """
        
        self._addProp("base_dir",base_dir,blob.FileBlob)
        self._readJson(self,os.path.join(self.getProperty("base_dir"),"info.json"))

 
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

        for k in json_in:
            self._addProp(k,json_in[k])
    
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
            

    @property
    def data(self):
        """
        Access this dataset in a programatically accessible way.
        """

        return None
        
    @property 
    def summary(self):
        """
        Return string summary of what this processing step actually did to the
        data.
        """

        output = ["Instance of '{:s}' class with the properties:".format(self.__name__)]
        for p in self._props:
            output.append("    {:s} : {}".format(p,self._props[p].value))

        
        return "\n".join(output) 
