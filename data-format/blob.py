__description__ = \
"""
Generic Blob classes to hold onto key/value pairs with specific validation.
"""
__author__ = "Michael J. Harms"
__date__ "2015-05-16"


class Blob:
    """
    Structure to hold arbitrary blobs of data.
    """

    def __init__(self,key_name,value):
        """
        Key and value.
        """
        self.key = key_name
        self.value = value
       
    def _validate_key(self):
        """
        Key validation.  Throw exception if bad.
        """
        pass    
     
    def _validate_value(self):
        """
        Value validation. Throw exception if bad.
        """
        pass        

    @property
    def key(self):
        """
        Return the key.
        """
        return self._key

    @key.setter
    def key(self,key):
        """
        Set the key after validation.
        """
        self._validate_key(key)
        self._key = key

    @property:
    def value(self):
        """
        Return the value.
        """
        return self._value

    @value.setter
    def value(self,value):
        """
        Set the value after validation. 
        """
        self._validate_value(value)
        self._value = value

class FileBlob:
    """
    Blob that checks to make sure the file added as a value is real.
    """
    def _validate_value(self,value):
        """
        Validate the file is there.  Throw exception if it can't open.
        """
   
        if not os.path.exists(value):
            err = "{:s} does not exist.\n".format(value)
            raise ValueError(err)
 
class DateBlob:
    """
    Blob that makes sure the date format in the value is valid.
    """
    
    def _validate_value(self,value):
        """
        Make sure the date is intelligible.
        """

        try: 
            datetime.datetime.strptime(value, "%Y-%m-%d")
        except (ValueError,TypeError):
            err = "Property for \"{}\" \"{}\" invalid.".format(key,value)
            raise ValueError(err)


