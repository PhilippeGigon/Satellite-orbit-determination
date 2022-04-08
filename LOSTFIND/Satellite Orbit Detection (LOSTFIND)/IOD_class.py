import math
import os
import re
import numpy as np
import datetime
import constants
from Coordinates_conversion import *
#############################################################
# The class IOD encodes the IODs which have to be placed in
# the same folder as the script. The function Create_IODs reads
# the txt file and returns three instances of IODs corresponding
# to the tree lines in the txt file. Note that the unknown values of
# IOD parameters must be replaced by a placeholder, i.g. by 0, s.t. the
# form of the IOD is always the same.
#############################################################


def is_IOD(IOD_candidate):
    ''''Tests if a given string is an IOD with the right format or not'''
    match = re.match(
        "[0-9]{5}\s[0-9]{2}\s[0-9]{3}[A-Z]{3}\s[0-9]{4}\s[A-Z]\s[0-9]{17}\s[0-9]{2}\s[0-9]{2}\s[0-9]{7}[+-][0-9]{6}\s[0-9]{2}\s[A-Z]\+[A-Z a-z]{3}\s[A-Z a-z]{2}\s[A-Z a-z]{6}", IOD_candidate)
    return match


def Create_IODs():
    ''' Opens the txt file and creates three instances of
     IOD classes containing the IODs'''

    script_dir = os.path.dirname(__file__)  # Location of python script
    rel_path = "IOD_ELFIND2.txt"  # Name of .txt
    abs_file_path = os.path.join(script_dir, rel_path)
    Data = open(abs_file_path, "r")
    # Checks that we have 3 IOD's in the file
    line_count = 0
    for line in Data:
        if line != "\n":
            line_count += 1
    if line_count != 3:
        Data.close()
        raise ValueError("There must be three IOD's in ", rel_path)
    if line_count == 3:
        Data.seek(0)
        IODs = Data.readlines()
        IODset = []
        for i in IODs:
            if is_IOD(i):
                IODset.append(IOD(i))
        Data.close()
        if len(IODset) != 3:
            raise ValueError("IODs have wrong format")
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

        # Computes UNIX epoch time
        self.Date = f[0:8]
        # Time
        t_year = int(f[0:4])
        t_month = int(f[4:6])
        t_day = int(f[6:8])
        t_hour = int(f[8:10])
        t_minute = int(f[10:12])
        t_second = int(f[12:14])
        t_microsecond = int(f[14:17])*1000  # in micro-second
        # in miliseconds
        # time in seconds since 01.01.1970 (UNIX)
        self.unix_epochtime = datetime.datetime(
            t_year, t_month, t_day, t_hour, t_minute, t_second, t_microsecond).timestamp()
        ###############################################
        ###############################################

        # Measured angles
        # Right ascention declination
        if (int(h[0:1])==1 or int(h[0:1])==2 or int(h[0:1])==3 or int(h[0:1])==7):
            self.type = 1
        else:
            self.type = 0

        if (int(h[0:1])==1):
            X, Y = i.split("-")
            RaH = int(X[0:2])
            RaM = int(X[2:4])
            RaS = int(X[4:6])
            Ras = int(X[6:7])
            self.RaRad = ((360/24)*RaH+15*RaM/60+15*RaS/3600 +
                          15*Ras/36000)*constants.degtorad  # Radians
            DecD = int(Y[0:2])
            DecM = int(Y[2:4])
            DecS = int(Y[4:6])
            if "-" in i:
                 DecRad = -(DecD + DecM/60+DecS/3600) * \
                     constants.degtorad  # Radians
            else: 
                self.DecRad = (DecD + DecM/60+DecS/3600) * \
                     constants.degtorad  # Radians

        if (int(h[0:1])==2):
            X, Y = i.split("-")
            RaH = int(X[0:2])
            RaM = int(X[2:4])
            Ram = int(X[4:7])
            self.RaRad = ((360/24)*RaH+15*RaM/60+15*Ram/60000)*constants.degtorad  # Radians
            DecD = int(Y[0:2])
            DecM = int(Y[2:4])
            Decm = int(Y[4:6])
            if "-" in i:
                self.DecRad = -(DecD + DecM/60+Decm/6000) * \
                     constants.degtorad  # Radians
            else:
                self.DecRad = (DecD + DecM/60+Decm/6000) * \
                     constants.degtorad  # Radians


        if (int(h[0:1])==3):
            X, Y = i.split("-")
            RaH = int(X[0:2])
            RaM = int(X[2:4])
            Ram = int(X[4:7])
            self.RaRad = ((360/24)*RaH+15*RaM/60+15*Ram/60000)*constants.degtorad  # Radians
            DecD = int(Y[0:2])
            Decd = int(Y[2:6])
            if "-" in i:
                self.DecRad = -(DecD + Decd/1000) * \
                     constants.degtorad  # Radians
            else:
                self.DecRad = (DecD + Decd/1000) * \
                     constants.degtorad  # Radians

        if (int(h[0:1])==7):
            X, Y = i.split("-")
            RaH = int(X[0:2])
            RaM = int(X[2:4])
            RaS = int(X[4:6])
            Ras = int(X[6:7])
            self.RaRad = ((360/24)*RaH+15*RaM/60+15*RaS/3600 +
                          15*Ras/36000)*constants.degtorad  # Radians
            DecD = int(Y[0:2])
            Decd = int(Y[2:6])
            if "-" in i:
                 self.DecRad = -(DecD + Decd/1000) * \
                         constants.degtorad  # Radians
            else:
                 self.DecRad = (DecD + Decd/1000) * \
                         constants.degtorad  # Radians

        if (int(h[0:1])==4):
            X, Y = i.split("-")
            AzD = int(X[0:3])
            AzM = int(X[3:5])
            AzS = int(X[5:7])
            self.AzRad = (AzD+AzM/60+AzS/3600)*constants.degtorad  # Radians
            ElD = int(Y[0:2])
            ElM = int(Y[2:4])
            ElS = int(Y[4:6])
            if "-" in i:
                 self.ElRad = -(ElD + ElM/60+ElS/3600) * \
                         constants.degtorad  # Radians
            else:
                 self.ElRad = (ElD + ElM/60+ElS/3600) * \
                         constants.degtorad  # Radians

        if (int(h[0:1])==5):
            X, Y = i.split("-")
            AzD = int(X[0:3])
            AzM = int(X[3:5])
            Azm = int(X[5:7])
            self.AzRad = (AzD+AzM/60+Azm/6000)*constants.degtorad  # Radians
            ElD = int(Y[0:2])
            ElM = int(Y[2:4])
            Elm = int(Y[4:6])
            if "-" in i:
                 self.ElRad = -(ElD + ElM/60+Elm/6000) * \
                         constants.degtorad  # Radians
            else:
                 self.ElRad = (ElD + ElM/60+Elm/6000) * \
                         constants.degtorad  # Radians

        if (int(h[0:1])==6):
            X, Y = i.split("-")
            AzD = int(X[0:3])
            Azd = int(X[3:7])
            self.AzRad = (AzD+Azd/10000)*constants.degtorad  # Radians
            ElD = int(Y[0:2])
            Eld = int(Y[2:6])
            if "-" in i:
                 self.ElRad = -(ElD + Eld/1000) * \
                         constants.degtorad  # Radians
            else:
                 self.ElRad = (ElD +Eld/10000) * \
                         constants.degtorad  # Radians
        #if "-" in i:
        #    X, Y = i.split("-")
        #    RaH = int(X[0:2])
        #    RaM = int(X[2:4])
        #   RaS = int(X[4:6])
        #  Ras = int(X[6:7])
        #    self.RaRad = ((360/24)*RaH+15*RaM/60+15*RaS/3600 +
        #                  15*Ras/36000)*constants.degtorad  # Radians
        #    DecD = int(Y[0:2])
        #    DecM = int(Y[2:4])
        #    DecS = int(Y[4:6])
        #    self.DecRad = -(DecD + DecM/60+DecS/3600) * \
        #        constants.degtorad  # Radians
            #print(self.RaRad/constants.degtorad)
            #print(self.DecRad/constants.degtorad)

        #if "+" in i:
        #    X, Y = i.split("+")
        #    RaH = int(X[0:2])
        #    RaM = int(X[2:4])
        #    RaS = int(X[4:6])
        #    Ras = int(X[6:7])
        #    self.RaRad = ((360/24)*RaH+RaM/60+RaS/3600 +
        #                  Ras/36000)*constants.degtorad  # Radians
        #    DecD = int(Y[0:2])
        #    DecM = int(Y[2:4])
        #    DecS = int(Y[4:6])
        #    self.DecRad = (DecD + DecM/60+DecS/3600) * \
        #        constants.degtorad  # Radians

    def get_e(self, latitude, time):
        '''returns the unit vector pointing from telescope towards satellite'''
        if(self.type == 1):
            alt = get_alt(self.DecRad,self.RaRad,latitude,time)
            az = get_az(alt,self.DecRad,self.RaRad,latitude,time)
        else:
            alt = self.ElRad
            az = self.AzRad
        e = np.array([math.cos(alt)*math.sin(az), math.cos(alt)*math.cos(az), math.sin(alt)]) 
        Rot = np.array([[-math.sin(time), -math.sin(latitude)*math.cos(time), math.cos(latitude)*math.cos(time)], [math.cos(time), -math.sin(latitude)*math.sin(time), math.cos(latitude)*math.sin(time)], [0, math.cos(latitude), math.sin(latitude)]])
        eloc = Rot.dot(e)
        return eloc

    def get_time(self):
        return self.unix_epochtime
 
    #def assign_az_el(ra_g,dec_g):
    #    self.ElRad = el;
    #    self.AzRad = az;
 

    def printIOD(self):
        print("Ra: ", self.RaRad, "Dec: ", self.DecRad)
        print("TIME SINCE UNIX:", self.unix_epochtime)


#############################################################
