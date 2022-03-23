from numpy import *
import os
import ephem
import datetime
import math
from IOD_class import *
from StationPositionVector import *
from UniversalKeplerEquation import *
from constants import mu


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
    
    # Extracting observation time and time intervalls:
    t1 = iodset[0].get_time()
    t2 = iodset[1].get_time()
    t3 = iodset[2].get_time()
    time_s1 = get_lst(lon1,t1) 
    time_s2 = get_lst(lon2,t2) 
    time_s3 = get_lst(lon3,t3) 
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

    ##########################################
    #################TEST#####################
    ##########################################
    R1 = np.array([3489800, 3430200, 4078500])
    R2 = np.array([3460100, 3460100, 4078500])
    R3 = np.array([3429900, 3490100,  4078500])
    e1 = np.array([0.71643, 0.68074, -0.15270])
    e2 = np.array([0.56897, 0.79531, -0.20917])
    e3 = np.array([0.41841, 0.87007, -0.26059])
    tau = 237.58
    tau1 = -118.1
    tau3 = 119.47
    ##########################################
    #################TEST#####################
    ##########################################

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
    b = -2*mu*B*(A+E)
    c = -mu**2*B**2

    # Norm of r2
    r2_norm = solveroot(a, b, c)

    rho1 = 1/D*((6*(D31*tau1/tau3+D21*tau/tau3)*r2_norm**3+mu*D31 *
                 (tau**2-tau1**2)*tau1/tau3)/(6*r2_norm**3+mu*(tau**2-tau3**2))-D11)
    rho2 = A+mu*B/r2_norm**3

    rho3 = 1/D*((6*(D13*tau3/tau1-D23*tau/tau1)*r2_norm**3+mu*D13 *
                 (tau**2-tau3**2)*tau3/tau1)/(6*r2_norm**3+mu*(tau**2-tau1**2))-D33)

    # Vectors pointing to satellite
    r1 = rho1*e1+R1
    r2 = rho2*e2+R2
    r3 = rho3*e3+R3

    # Lagrange coeffs
    f1 = 1-1/2*mu/r2_norm**3*tau1**2
    g1 = tau1-1/6*mu/r2_norm**3*tau1**3
    f3 = 1-1/2*mu/r2_norm**3*tau3**2
    g3 = tau3-1/6*mu/r2_norm**3*tau3**3

    # speed at position 2
    v2 = 1/(f1*g3-f3*g1)*(-f3*r1+f1*r3)

    # iterate:
    r2 = np.array([5659100, 6533800, 3270100])
    v2 = np.array([-3908, 5057.3, -2222.2])
    chi1 = get_chi(tau1, r2/1000, v2/1000)
    chi3 = get_chi(tau3, r2/1000, v2/1000)
    print(chi1, chi3)
    return r1, r3


IODset = array(Create_IODs())
r1, r2 = find_r(IODset)
