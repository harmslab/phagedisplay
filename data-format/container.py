
class Container:
    """
    """

    def __init__(self):
        """
        """
   
        # ._data will hold instnaces of sub-classes of ProcessorBase
        self._data = []

    def addStep(self,process_class):

        self._data.append(process_class.process(self.data))

    def _write(self):
        """
        """

        pass   

    @property
    def current_data(self):
        """
        """
        
        if len(self._data) == 0:
            return None
        
        return self._data[-1].data

