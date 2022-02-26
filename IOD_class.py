import math
import os

#############################################################
# The class IOD encodes the IODs which have to be placed in
# the same folder as the script. The function Create_IODs reads
# the txt file and returns three instances of IODs corresponding
# to the tree lines in the txt file. Note that the unknown values of
# IOD parameters must be replaced by a placeholder, i.g. by 0, s.t. the
# form of the IOD is always the same.
#############################################################


def Create_IODs():
    ''' Opens the txt file and creates three instances of
     IOD classes containing the IODs'''

    script_dir = os.path.dirname(__file__)  # Location of python script
    rel_path = "test.txt"  # Name of .txt
    abs_file_path = os.path.join(script_dir, rel_path)

    Data = open(abs_file_path, "r")
    IODs = Data.readlines()
    IODset = []
    for i in IODs:
        IODset.append(IOD(i))
    Data.close()
    return IODset

#############################################################


class IOD:
    # Encodes the initial data (IOD)
    def __init__(self, IOD_string):
        # The different elements are seperated by empty space in IOD
        # Methods
        a, b, c, d, e, f, g, h, i, j, k, l, m = IOD_string.split(" ")
        self.Object_ID = a
        self.Station_ID = b
        self.Station_stat = c
        # Date year month day
        self.Date = f[0:8]
        # Time
        self.t_hour = int(f[8:10])
        self.t_minute = int(f[10:12])
        self.t_msecond = int(f[12:])  # in miliseconds
        # Measured angles
        # El either positive or negative
        if "-" in i:
            X, Y = i.split("-")
            self.Az = int(X)
            self.El = -int(Y)
        if "+" in i:
            X, Y = i.split("+")
            self.Az = int(X)
            self.El = int(Y)

    def get_Az(self):
        return self.Az

    def get_El(self):
        return self.El

    def get_miliseconds(self):
        return self.t_msecond

    def printIOD(self):
        print("Az: ", self.Az, "El: ", self.El)
        print("HH: ", self.t_hour, "MM: ",
              self.t_minute, "SSSS: ", self.t_msecond)

#############################################################


# Debbug
iodset = Create_IODs()
for i in iodset:
    i.printIOD()
