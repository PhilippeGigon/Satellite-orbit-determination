import math
import constants
import numpy as np

class TLE: 
    def __init__(self, input, precision):
        self.ra = np.array(input[0])
        self.rb = input[1]
        self.ta = input[2]
        self.tb = input[3]
        self.prec = precision
        self.W = self.unit_vect_W()

        # calculated TLE components
        self.e = self.eccentricity()
        self.n = self.mean_motion()
        self.i = self.inclination(self.W)
        self.Omega = self.right_ascension(self.W)
        self.omega = self.arg_perigee()
        self.M = self.mean_anomaly()
        self.t = (self.ta + self.tb)/2

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





    def __getitem__(self, index):
        return self.W[index]
    
    
    #calculates the angle between the two vectors
    def true_anomaly(self):
        n = np.linalg.norm(self.ra)*np.linalg.norm(self.rb)
        d = np.dot(self.ra, self.rb)
        return np.arccos(d/n)

    def rad2deg(self, arg):
        argDeg = arg*360/(2*constants.pi)
        if argDeg < 0:
            argDeg += 360
        return argDeg

    #calculates the triangle area
    def calc_delta(self): 
        nra = np.linalg.norm(self.ra)
        nrb = np.linalg.norm(self.rb)
        drab = np.dot(self.ra,self.rb)
        # triangle delta
        dnu = self.true_anomaly()
        delta = 0.5 * nra * nrb * np.sin(dnu)
    
        return delta


    # Two functions W and f used in calc_eta
    def fct_W(self, w):
        g = 2*np.arcsin(np.sqrt(w))
        return (2*g - np.sin(2*g))/(np.sin(g)**3)

    def fct_f(self, x, m, l):
        W = self.fct_W(m/x**2 - l)
        return 1 - x + m/x**2 * W 


    # Calculates eta = S/delta
    def calc_eta(self):
        nra = np.linalg.norm(self.ra)
        nrb = np.linalg.norm(self.rb)
        drab = np.dot(self.ra,self.rb)

        tau = np.sqrt(constants.G * constants.MT) * (self.tb - self.ta)
        m = tau**2 / (2*nra*nrb + drab)**(3/2)
        l = (nra + nrb)/(2*np.sqrt(2*(nra*nrb + drab))) - 0.5
    
        # iteself.self.rative process to find eta
        eta0 = 12/22 + 10/22*np.sqrt(1+44/9*m/(l+(5/6)))
        eta1 = eta0 + 0.1
        eta2 = eta0
        f = self.fct_f(eta2, m, l)
        while f > self.prec:
            eta2 = eta1 - self.fct_f(eta1,m,l) * (eta1-eta0) / (self.fct_f(eta1,m,l)-self.fct_f(eta0,m,l))
            f = self.fct_f(eta2,m,l)
            eta0 = eta1
            eta1 = eta2
        eta = eta2
        return eta

    # unit vector ea x e0
    def unit_vect_W(self):
        ea = (1/np.linalg.norm(self.ra))*self.ra
        r0 = self.rb - (np.dot(self.rb,ea))*ea
        e0 = (1/np.linalg.norm(r0))*r0
        W = np.cross(ea, e0)
        return W

    def calc_nua(self):
        delta = self.calc_delta()
        eta = self.calc_eta()
        dnu = self.true_anomaly()

        tau = np.sqrt(constants.G * constants.MT) * (self.tb - self.ta)
        p = (2*delta*eta/tau)**2
        pa = p/np.linalg.norm(self.ra) -1
        pb = p/np.linalg.norm(self.rb) -1
        nua = np.arctan((pa*np.cos(dnu) - pb))
        return nua


    ##################### 
    # Functions to calculate the TLE components
    #####################

    
    def inclination(self, W):
        arg = np.sqrt(W[0]**2 + W[1]**2)/W[2]
        inc = np.arctan(arg)
        return self.rad2deg(inc)

    def right_ascension(self, W):
        Omega = np.arctan(W[0]/-W[1])
        return self.rad2deg(Omega)

    def arg_perigee(self):
        W = self.unit_vect_W()
        nua = self.calc_nua()
        ua = np.arctan(self.ra[2]/(-self.ra[0]*W[1] + self.ra[1]*W[1]))
        arg = ua - nua
        return self.rad2deg(arg)

    def eccentricity(self):
        delta = self.calc_delta()
        eta = self.calc_eta()
        dnu = self.true_anomaly()

        tau = np.sqrt(constants.G * constants.MT) * (self.tb - self.ta)
        p = (2*delta*eta/tau)**2
        pa = p/np.linalg.norm(self.ra) -1
        pb = p/np.linalg.norm(self.rb) -1
        nua = np.arctan((pa*np.cos(dnu) - pb)/(pa*np.sin(dnu)))
        e = pa/np.cos(nua)
        return e

    def mean_motion(self):
        delta = self.calc_delta()
        eta = self.calc_eta()
        tau = np.sqrt(constants.G * constants.MT) * (self.tb - self.ta)
        p = (2*delta*eta/tau)**2
        a = p/(1-self.e**2)
        n = np.sqrt(constants.G*constants.MT/(a**3))
        return n

    def mean_anomaly(self):
        e = self.eccentricity()
        nua = self.calc_nua()
        Ea = np.arctan((np.sqrt(1-e**2)*np.sin(nua)) / (np.cos(nua) + e))
        M = Ea- e*np.sin(Ea)
        return self.rad2deg(M)

    #################################
    # Additional TLE components
    #################################

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


    # create list of TLE
    def list_TLE(self):
        return [self.e, self.n, self.i, self.Omega, self.omega, self.M, self.t]

    def TLE_format(self):
        line1 = str(1) + ' ' + self.SCN + self.Class + ' ' + self.intDesY \
            + self.intDesNY + self.intDesP + ' ' + self.epochY + self.epochD \
            + ' ' + self.derMM + ' ' + self.der2MM + ' ' + self.Bstar + ' ' \
            + self.EphType + ' ' + self.ESN + str(self.checksum)
        # inclination
        str_i = f"{self.i:.4f}"
        if self.i < 100:
            str_i = '0' + str_i
    
        # right ascension
        str_O = "{:8.4F}".format(self.Omega)
        # eccentricity
        string_e = str(round(self.e, 7))
        str_e = string_e[3:]
        # argument of perigee
        str_o = "{:8.4F}".format(self.omega)
        # mean anomaly
        str_M = str(round(self.M, 4))
        if self.M < 100:
            str_M = '0' + str_M
        # mean motion
        str_n = str(round(self.n, 8))
        if self.n < 10:
            str_n = '0' + str_n
        # revolution 
        str_rev = "{:5d}".format(self.rev)
        

        line2 = str(2) + ' ' + self.SCN + ' ' + str_i + ' ' + str_O + ' '\
            + str_e + ' ' + str_o + ' ' + str_M + ' ' + str_n + str_rev \
            + str(self.checksum)

        f = open("TLE.out", "w")
        f.write(line1 + "\n" + line2)
        f.close