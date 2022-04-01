import math
import numpy as np
import constants

class TLE:
    def __init__(self, input):
        self.r = input[0]
        self.v = input[1]
        self.t = input[2]

         # additional TLE components
        self.SCN = self.getSCN() # Satellite Catalog Number
        self.Class = self.getClass() # classification of the satellite U = unclassified
        self.intDesY = self.getIntDesY() # International Designator = last two digits of launch year
        self.intDesNY = self.getIntDesNY() # International Designator = launch number of the year
        self.intDesP = self.getIntDesP() # International Designator = piece of the launch
        self.epochY = self.getEpochY() # last two digits of epoch year
        self.epochD = self.getEpochD() # day and fraction of the day
        self.derMM = self.getDerMM() # first derivative of the mean motion
        self.der2MM = self.getDer2MM() # second derivative of the mean motion
        self.Bstar = self.getBstar() # Drag term, radiation pressure coeff
        self.EphType = self.getEphType() # Ephemeris type, always 0
        self.ESN = self.getESN() # Element set number
        self.checksum = self.getChecksum() # Checksum
        self.rev = self.getRev() # revolution number at epoch



    def values(self):
        h = np.cross(self.r, self.v)
        z = [0,0,1]
        N = np.cross(z,h)
        rn = np.linalg.norm(self.r)
        vr = np.dot(self.r, self.v)/rn

        # inclination
        i = np.arccos(h[2]/np.linalg.norm(h))

        # right ascension
        Omega = np.arccos(N[2]/np.linalg.norm(N))

        #eccentricity
        ev = 1/constants.mukm * ((rn**2 - constants.mukm/rn)*self.r - rn*vr*self.v)
        e = np.linalg.norm(ev)
        
        # argument of perigee
        omega = np.arccos(np.dot(N, ev)/(e*np.linalg.norm(N)))

        # true anomaly 
        theta = np.arccos(np.dot(ev, self.r)/(e*rn))

        # mean anomaly
        E = np.arctan((np.sqrt(1-e**2)*np.sin(theta)) / (np.cos(theta) + e))
        M = E - e*np.sin(E)

        # mean motion
        a = np.linalg.norm(h)**2/(constants.mu*(1-e**2))
        n = np.sqrt(constants.mu/a**3)

        return i, Omega, e, omega, M, n

    def getSCN(self):
        SCN = 0
        return f"{SCN:05d}"
    
    def getClass(self):
        Class = 'U'
        return Class

    def getIntDesY(self):
        IntDesY = 0
        return f"{IntDesY:02d}"

    def getIntDesNY(self):
        IntDesNY = 0
        return f"{IntDesNY:03d}"

    def getIntDesP(self):
        IntDesP = 0
        return f"{IntDesP:03d}"
    
    def getEpochY(self):
        time = self.t
        years = int(time / (60*60*24*365)) + 1970
        if years < 2000:
            EpochY = years - 1900
        else: 
            EpochY = years - 2000
        return f"{EpochY:02d}"

    def getEpochD(self):
        time = self.t
        years = int(time / (60*60*24*365))
        EpochD = (time - (years * 60*60*24*365)) / (60*60*24)
        retprov = "{:11.8F}".format(EpochD)
        return retprov

    def getDerMM(self):
        DerMM = 0
        return "{:.7f}".format(DerMM)

    def getDer2MM(self):
        Der2MM = 0
        return f"{Der2MM:010d}"

    def getBstar(self):
        inp = 0
        #exp = fexp(inp)
        #first = int(inp * 10**(-exp+4))
        #second = exp + 1
        #return (str(first) + str(second))
        return f"{inp:08d}"

    def getEphType(self):
        EphType = '0'
        return EphType
    
    def getESN(self):
        ESN = 0
        return f"{ESN:02d}"

    def getChecksum(self):
        orbnum = 0
        return (orbnum % 10)

    def getRev(self):
        rev = 0
        return 0


    def TLE_format(self):
        i, Omega, e, omega, M, n = self.values()

        line1 = str(1) + ' ' + self.SCN + self.Class + ' ' + self.intDesY \
            + self.intDesNY + self.intDesP + ' ' + self.epochY + self.epochD \
            + ' ' + self.derMM + ' ' + self.der2MM + ' ' + self.Bstar + ' ' \
            + self.EphType + ' ' + self.ESN + str(self.checksum)
        
        # inclination
        str_i = f"{i:.4f}"
        if i < 100:
            str_i = '0' + str_i
        # right ascension
        str_O = "{:8.4F}".format(Omega)
        # eccentricity
        string_e = str(round(e, 7))
        str_e = string_e[3:]
        # argument of perigee
        str_o = "{:8.4F}".format(omega)
        # mean anomaly
        #str_M = str(round(M, 4))
        str_M = "{:11.8F}".format(M)
        #if M < 100:
        #    str_M = '0' + str_M
        # mean motion
        str_n = str(round(n, 8))
        if n < 10:
            str_n = '0' + str_n
        # revolution 
        str_rev = "{:5d}".format(self.rev)

        line2 = str(2) + ' ' + self.SCN + ' ' + str_i + ' ' + str_O + ' '\
            + str_e + ' ' + str_o + ' ' + str_M + ' ' + str_n + str_rev \
            + str(self.checksum)

        f = open("TLE.out", "w")
        f.write(line1 + "\n" + line1)
        f.close