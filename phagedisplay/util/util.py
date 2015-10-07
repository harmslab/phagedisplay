
import datetime

def logger(message,output_file):
    """
    Log a message to a file, appending the date.
    """

    d = "_".join(str(datetime.date.today()).split(" "))

    f = open(output_file,'a')
    f.write(format("{:25s}{:s}\n",d,message))
    f.close()
