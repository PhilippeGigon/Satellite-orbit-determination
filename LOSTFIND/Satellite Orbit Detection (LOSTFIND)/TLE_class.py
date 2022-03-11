import math
import constants
import numpy as np

class TLE: 
    def __init__(self, input, precision):
        self.self.self.ra = input[0]
        self.rb = input[1]
        self.ta = input[2]
        self.tb = input[3]
        self.prec = precision
        self.e = self.eccentricity(self)
        self.n = 0 ##ToDo
        self.i = self.inclination(self)
        self.Omega = self.right_ascension(self)
        self.omega = self.arg_perigee(self)
        self.M = self.mean_anomaly(self)
        self.t = 0 ##ToDo
    
    
    #calculates the angle between the two vectors
    def true_anomaly(self):
        n = np.linalg.norm(self.ra)*np.linalg.norm(self.rb)
        d = np.linalg.dot(self.ra, self.rb)
        return np.arccos(d/n)


    #calculates the triangle area
    def calc_delta(self): 
        nra = np.linalg.norm(self.ra)
        nrb = np.linalg.norm(self.rb)
        drab = np.linalg.dot(self.ra,self.rb)
        # triangle delta
        dnu = self.true_anomaly(self.ra, self.rb)
        delta = 0.5 * nra * nrb * np.sin(dnu)
    
        return delta


    # Two functions W and f used in calc_eta
    def fct_W(w):
        g = 2*np.arcsin(np.sqrt(w))
        return (2*g - np.sin(2*g))/(np.sin(g)**3)

    def fct_f(self, x, m, l):
        W = self.fct_W(m/x**2 - l)
        return 1 - x + m/x**2 * W 


    # Calculates eta = S/delta
    def calc_eta(self):
        nra = np.linalg.norm(self.ra)
        nrb = np.linalg.norm(self.rb)
        drab = np.linalg.dot(self.ra,self.rb)

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

    def unit_vect_W(self):
        ea = (1/np.linalg.norm(self.ra))*self.ra
        r0 = self.rb - (np.linalg.dot(self.rb,ea))*ea
        e0 = (1/np.linalg.norm(r0))*r0
        W = np.linalg.cross(ea, e0)
        return W

    def calc_nua(self):
        delta = self.calc_delta(self)
        eta = self.calc_eta(self)
        dnu = self.true_anomaly(self)

        tau = np.sqrt(constants.G * constants.MT) * (self.tb - self.ta)
        p = (2*delta*eta/tau)**2
        pa = p/np.linalg.norm(self.ra) -1
        pb = p/np.linalg.norm(self.rb) -1
        nua = np.arctan((pa*np.cos(dnu) - pb))
        return nua


    ##################### 
    # Functions to calculate the TLE components
    #####################

    
    def inclination(W):
        arg = np.sqrt(W[0]**2 + W[1]**2)/W[2]
        return np.arctan(arg)

    def right_ascension(W):
        return np.arctan(W[0]/-W[1])

    def arg_perigee(self):
        W = self.unit_vect_W(self)
        nua = self.calc_nua(self)
        ua = np.arctan(self.ra[2]/(-self.ra[0]*W[1] + self.ra[1]*W[1]))
        return ua - nua

    def eccentricity(self):
        delta = self.calc_delta(self)
        eta = self.calc_eta(self)
        dnu = self.true_anomaly(self)

        tau = np.sqrt(constants.G * constants.MT) * (self.tb - self.ta)
        p = (2*delta*eta/tau)**2
        pa = p/np.linalg.norm(self.ra) -1
        pb = p/np.linalg.norm(self.rb) -1
        nua = np.arctan((pa*np.cos(dnu) - pb))
        e = pa/np.cos(nua)
        return e

    def mean_anomaly(self):
        e = self.eccentricity(self)
        nua = self.calc_nua(self)
        Ea = np.arctan((np.sqrt(1-e**2)*np.sin(nua)) / (np.cos(nua) + e))
        M = Ea- e*np.sin(Ea)
        return M

    # create list of TLE
    def list_TLE(self):
        return [self.e, self.n, self.i, self.Omega, self.omega, self.M, self.t]
