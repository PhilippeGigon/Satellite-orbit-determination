import math
import constants
import numpy as np
from delta_eta import *

# Calculates the unit vectors ea and e0 and returns their cross product W
def unit_vect_W(ra,rb):
    ea = (1/np.linalg.norm(ra))*ra
    r0 = rb - (np.linalg.dot(rb,ea))*ea
    e0 = (1/np.linalg.norm(r0))*r0
    W = np.linalg.cross(ea, e0)
    return W

def inclination(W):
    arg = np.sqrt(W[0]**2 + W[1]**2)/W[2]
    return np.arctan(arg)

def right_ascension(W):
    return np.arctan(W[0]/-W[1])

def arg_perigee(ra, nua, W):
    ua = np.arctan(ra[2]/(-ra[0]*W[1] + ra[1]*W[1]))
    return ua - nua

def eccentricity(ra, rb,  ta, tb, prec):
    delta = calc_delta(ra, rb, ta, tb)
    eta = calc_eta(ra, rb, ta, tb, prec)
    dnu = true_anomaly(ra, rb)

    tau = np.sqrt(constants.G * constants.MT) * (tb - ta)
    p = (2*delta*eta/tau)**2
    pa = p/np.linalg.norm(ra) -1
    pb = p/np.linalg.norm(rb) -1
    nua = np.arctan((pa*np.cos(dnu) - pb))
    e = pa/np.cos(nua)
    return e


