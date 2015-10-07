
import datetime, sys

def logger(message,output_file):
    """
    Log a message to a file, appending the date.
    """

    d = "_".join(str(datetime.date.today()).split(" "))

    f = open(output_file,'a')
    f.write("{:25s}{:s}\n".format(d,message))
    f.close()

    print(message)
    sys.stdout.flush()
