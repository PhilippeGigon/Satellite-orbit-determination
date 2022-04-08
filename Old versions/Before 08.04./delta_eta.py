import math
import constants
import numpy as np 



#calculates the angle between the two vectors
def true_anomaly(ra, rb):
    n = np.linalg.norm(ra)*np.linalg.norm(rb)
    d = np.linalg.dot(ra, rb)
    return np.arccos(d/n)


#calculates the triangle area
def calc_delta(ra, rb, ta, tb): 
    nra = np.linalg.norm(ra)
    nrb = np.linalg.norm(rb)
    drab = np.linalg.dot(ra,ra)
    # triangle delta
    dnu = true_anomaly(ra, rb)
    delta = 0.5 * nra * nrb * np.sin(dnu)
    
    return delta


# Two functions W and f used in calc_eta
def fct_W(w):
    g = 2*np.arcsin(np.sqrt(w))
    return (2*g - np.sin(2*g))/(np.sin(g)**3)

def fct_f(x, m, l):
    W = fct_W(m/x**2 - l)
    return 1 - x + m/x**2 * W 


# Calculates eta = S/delta
def calc_eta(ra, rb, ta, tb, prec):
    nra = np.linalg.norm(ra)
    nrb = np.linalg.norm(rb)
    drab = np.linalg.dot(ra,ra)

    tau = np.sqrt(constants.G * constants.MT) * (tb - ta)
    m = tau**2 / (2*(nra*nrb + drab))**(3/2)
    l = (nra + nrb)/(2*np.sqrt(2*(nra*nrb + drab))) - 0.5
    
    # iterative process to find eta
    eta0 = 12/22 + 10/22*np.sqrt(1+44/9*m/(l+(5/6)))
    eta1 = eta0 + 0.1
    eta2 = eta0
    f = fct_f(eta2, m, l)
    while f > prec:
        eta2 = eta1 - fct_f(eta1,m,l) * (eta1-eta0) / (fct_f(eta1,m,l)-fct_f(eta0,m,l))
        f = fct_f(eta2,m,l)
        eta0 = eta1
        eta1 = eta2
    eta = eta2
    return eta

