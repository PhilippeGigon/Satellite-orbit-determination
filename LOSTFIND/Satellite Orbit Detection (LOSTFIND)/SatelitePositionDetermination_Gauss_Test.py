from numpy import *
import os
import ephem
import datetime
import math
from IOD_class import *
from StationPositionVector import *
from UniversalKeplerEquation import *
from constants import mu, mukm


def solveroot(A, B, C):
    '''This function solves the root
    of the polynomial x^8+ax^6+bx^3+c=0 which
    corresponds to the norm of r2, it needs to keep only
    the physical solutions (positive and not complex)'''
    coefs = np.array([1, 0, A, 0, 0, B, 0, 0, C])
    root = np.roots(coefs)

    n_roots = 0
    for i in root:
        if (i > 0) and np.isreal(i):
            candidate = i.real
            n_roots = n_roots+1
    if n_roots == 0:
        raise ValueError("No root found")
    if n_roots > 1:
        raise ValueError("Multiple roots found")

    return candidate


def find_r(iodset):
    '''This code takes an array of IOD' and then computes
    two vectors r1 and r2 pointing to the satellite
    It uses the Gauss method'''

    # Extracting location of station
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

    t1 = iodset[0].get_time()
    t2 = iodset[1].get_time()
    t3 = iodset[2].get_time()
    time_s1 = get_lst(lon1, t1)
    time_s2 = get_lst(lon2, t2)
    time_s3 = get_lst(lon3, t3)

    ############################
    #############TEST###########
    ############################
    t1 = 0
    t2 = 118.10
    t3 = 237.58
    time_s1 = 44.506*constants.degtorad
    time_s2 = 45.0*constants.degtorad
    time_s3 = 45.499*constants.degtorad
    ############################
    #############TEST###########
    ############################
    tau = t3-t1
    tau1 = t1-t2
    tau3 = t3-t2

    # Creates the vector pointing to the station
    R1 = get_R(time_s1, lat1, h1)
    R2 = get_R(time_s2, lat2, h2)
    R3 = get_R(time_s3, lat3, h3)

    # Unit vectors pointing to satellite
    e1 = iodset[0].get_e(lat1, time_s1)
    e2 = iodset[1].get_e(lat2, time_s2)
    e3 = iodset[2].get_e(lat3, time_s3)

    # Usefull definitions
    p1 = cross(e2, e3)
    p2 = cross(e1, e3)
    p3 = cross(e1, e2)
    D = dot(e1, p1)

    D11 = dot(R1, p1)
    D12 = dot(R1, p2)
    D13 = dot(R1, p3)
    D21 = dot(R2, p1)
    D22 = dot(R2, p2)
    D23 = dot(R2, p3)
    D31 = dot(R3, p1)
    D32 = dot(R3, p2)
    D33 = dot(R3, p3)

    A = 1/D*(-D12*tau3/tau+D22+D32*tau1/tau)
    B = 1/(6*D)*(D12*(tau3**2-tau**2)*tau3/tau+D32*(tau**2-tau1**2)*tau1/tau)
    E = dot(R2, e2)
    R2normsquar = dot(R2, R2)

    a = -(A**2+2*A*E+R2normsquar)
    b = -2*mukm*B*(A+E)
    c = -mukm**2*B**2

    # Norm of r2
    r2_norm = solveroot(a, b, c)

    rho1 = 1/D*((6*(D31*tau1/tau3+D21*tau/tau3)*r2_norm**3+mukm*D31 *
                 (tau**2-tau1**2)*tau1/tau3)/(6*r2_norm**3+mukm*(tau**2-tau3**2))-D11)
    rho2 = A+mukm*B/r2_norm**3

    rho3 = 1/D*((6*(D13*tau3/tau1-D23*tau/tau1)*r2_norm**3+mukm*D13 *
                 (tau**2-tau3**2)*tau3/tau1)/(6*r2_norm**3+mukm*(tau**2-tau3**2))-D33)
    # Vectors pointing to satellite
    r1 = rho1*e1+R1
    r2 = rho2*e2+R2
    r3 = rho3*e3+R3
    print(r1, r2, r3)
    # Lagrange coeffs
    f1 = 1-1/2*mu/r2_norm**3*tau1**2
    g1 = tau1-1/6*mu/r2_norm**3*tau1**3
    f3 = 1-1/2*mu/r2_norm**3*tau3**2
    g3 = tau3-1/6*mu/r2_norm**3*tau3**3

    # speed at position 2
    v2 = 1/(f1*g3-f3*g1)*(-f3*r1+f1*r3)

    ##############################################################
    ########################Iterative part########################
    ##############################################################

    rho1old = 2*rho1
    # while abs(rho1old-rho1)>0.001:
    for i in range(0, 0):
        # All in km, km/s:
        # computes norm
        r_mag = math.sqrt(np.dot(r2, r2))
        v_mag = math.sqrt(np.dot(v2, v2))
        # Reciprocal of semi major axis
        alpha = 2/r_mag-v_mag**2/mukm
        chi1 = get_chi(tau1, r2, v2)
        chi3 = get_chi(tau3, r2, v2)
        f1 = 1-chi1**2/r_mag*C(alpha*chi1**2)
        g1 = tau1-1/sqrt(mukm)*chi1**3*S(alpha*chi1**2)
        f3 = 1-chi3**2/r_mag*C(alpha*chi3**2)
        g3 = tau3-1/sqrt(mukm)*chi3**3*S(alpha*chi3**2)

        c1 = g3/(f1*g3-f3*g1)
        c3 = -g1/(f1*g3-f3*g1)
        # Old values of rho for convergence test
        rho1old = rho1
        rho2old = rho2
        rho3old = rho3
        # Updated values of rhoi
        rho1 = 1/D*(-D11+1/c1*D21-c3/c1*D31)
        rho2 = 1/D*(-c1*D12+D22-c3*D32)
        rho3 = 1/D*(-c1/c3*D13+1/c3*D23-D33)
        r1 = rho1*e1+R1
        r2 = rho2*e2+R2
        r3 = rho3*e3+R3
        v2 = 1/(f1*g3-f3*g1)*(-f3*r1+f1*r3)
    # Converts back to SI units
    r1 = r1*1000
    r2 = r2*1000
    r3 = r3*1000
    v2 = v2*1000
    return r1, r3


IODset = array(Create_IODs())
r1, r2 = find_r(IODset)
