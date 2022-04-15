################################################################
# ElFind preliminary orbit determinatin algorithm, based on:
# Fundamentals of astrodynamics by Bate, Mueller and White
###############################################################
import math
import numpy
import numpy as np
from StationPositionVector import *
import constants
from IOD_class import *
###############################################################

def Ri_Li(iodset):
    '''Extracting location of station'''
    script_dir = os.path.dirname(__file__)  # Location of python script
    rel_path = "LOCATION.txt"  # Name of .txt
    abs_file_path = os.path.join(script_dir, rel_path)
    Data = open(abs_file_path, "r")
    Station_Coordinates = Data.readlines()
    Data.close()
    lat1 = float(Station_Coordinates[0])*(np.pi/180)
    lon1 = float(Station_Coordinates[1])*(np.pi/180)
    h1 = float(Station_Coordinates[2])
    lat2 = float(Station_Coordinates[3])*(np.pi/180)
    lon2 = float(Station_Coordinates[4])*(np.pi/180)
    h2 = float(Station_Coordinates[5])
    lat3 = float(Station_Coordinates[6])*(np.pi/180)
    lon3 = float(Station_Coordinates[7])*(np.pi/180)
    h3 = float(Station_Coordinates[8])

    # Creates the vector pointing to the station
    R1 = get_R(iodset[0].get_time(), lat1, h1)
    R2 = get_R(iodset[1].get_time(), lat2, h2)
    R3 = get_R(iodset[2].get_time(), lat3, h3)

    R = [R1, R2, R3]

    # Unit vectors pointing to satellite
    L1 = iodset[0].get_e(lat1, lon1, R1)
    L2 = iodset[1].get_e(lat2, lon2, R2)
    L3 = iodset[2].get_e(lat3, lon3, R3)

    L = [L1, L2, L3]
    return [R, L]

###############################################################


def procedure(iodset):

    # Extraction of Ri and Li
    RL = Ri_Li(iodset)
    R1 = Ri_Li[0][0]
    R2 = Ri_Li[0][1]
    R3 = Ri_Li[0][2]

    L1 = Ri_Li[1][0]
    L2 = Ri_Li[1][1]
    L3 = Ri_Li[1][2]

    # Extracting observation time:
    t1 = iodset[0].get_time()
    t2 = iodset[1].get_time()
    t3 = iodset[2].get_time()

    # Estimation of r2
    lowOrbit = np.linalg.norm(R2)/3
    r2 = R2 + lowOrbit*L2

    # Calculation of u
    mu = constants.mu
    u = mu/(np.linalg.norm(r2)**3)

    # calculation of v2
    v2 = [1, 1, 1]        # guess a first value
    fgv2 = fandg(r2, v2, t3-t2)
    v2 = 0

    # Calculation of f1, g1, f3, g3

    # Estimation of the magnitude of r2
    r2 = np.norm(R2) + 10

###############################################################

def fandg(rr1, vv1, delt, fr, gr):
    '''Calculates the f and g series'''
    r, u, p, q, p2, p4, q2, u2 = 0
    f, g = []

    r = np.linalg.norm(rr1)
    u = 1/(r**3)
    p = np.dot(rr1, vv1)/(r**2)
    q = np.dot(vv1, vv1)/(r**2) - u
    p2 = p**2
    p4 = p**4
    q2 = q**2
    u2 = u**2

    f[0] = 1
    f[1] = 0
    f[2] = -u/2
    f[3] = p*u/2
    f[4] = u*(u - 3*(5*p2-q))/24
    f[5] = -p*u*(u-7*p2+3*q)/8
    f[6] = -u*(u2 - 6*(35*p2-4*q)*u + 45*(21*p4-14*p2*q+q2))/720
    f[7] = p*u*(u2 - 2*(25*p2-7*q)*u + 5*(33*p4-30*p2*q+5*q2))/80
    f[8] = u*(u2*u - 9*(245*p2-13*q)*u2 + 27*(1925*p4-910*p2*q+41*q2)*u
              - 315*(429*p4*p2 - 495*p4*q + 135*p2*q2 - 5*q2*q))/40320

    g[0] = 0
    g[1] = 1
    g[2] = 0
    g[3] = -u/6
    g[4] = p*u/4
    g[5] = u*(u-9*(5*p2-q))/120
    g[6] = -p*u*(u-2*(7*p2-3*q))/24
    g[7] = -u*(u2-18*(35*p2-3*q)*u + 225*(21*p4-14*p2*q+q2))/5040
    g[8] = p*u*(u2-4*(25*p2-6*q)*u + 15*(33*p4-30*p2*q+5*q2))/320

    fr = 1 + delt*delt*(f[2] + delt*(f[3] + delt*(f[4] + delt*(f[5] +
                                                               delt*(f[6] + delt*(f[7] + delt*f[8]))))))

    gr = delt*(1+delt*delt*(g[3] + delt*(g[4] + delt*(g[5] +
                                                      delt*(g[6] + delt*(g[7] + delt*g[8]))))))

    return [fr, gr]

# Gaussian elimination
###############################################################

def rref(m, b):
    '''Gaussian elimination'''
    i, j, k, ROW = 0
    s = []
    bin, mult = 0

    # calculate scale factors
    for i in range(6):
        s[i] = abs(m[i][0])
        for j in range(1, 6):
            if s[i] < abs(m[i][j]):
                s[i] = abs(m[i][j])

    for j in range(5):
        ROW = j
        for i in range(j+1, 6):
            if abs(m[ROW][j])/s[ROW] < abs(m[i][j]/s[i]):
                ROW = i
        if ROW != j:
            for k in range(j, 7):   # swap rows
                bin = m[j][k]
                m[j][k] = m[ROW][k]
                m[ROW][k] = bin
            bin = s[j]  # swap scales
            s[j] = s[ROW]
            s[ROW] = bin
        # forward elimination

        for i in range(j+1, 6):
            mult = m[i][j]/m[j][j]
            for k in range(j+1, 7):
                m[i][k] = m[i][k] - mult * m[j][k]
            m[i][j] = 0

    # test for singular matrix
    bin = 1
    for i in range(6):
        bin = m[i][i]
    if bin == 0:
        print('Singular matix')
        exit(0)

    # back substitution
    b[5] = m[5][6]/m[5][5]
    for i in range(4, 0, -1):
        bin = 0
        for k in range(i+1, 6):
            bin = bin + m[i][k]*b[k]
        b[i] = (m[i][6] - bin)/m[i][i]
