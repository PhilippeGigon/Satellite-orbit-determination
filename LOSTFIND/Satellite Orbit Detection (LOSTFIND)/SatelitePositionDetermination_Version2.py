from IOD_class import *
from numpy import *
import os
import ephem
import datetime
import math
import numpy as np
from StationPositionVector import *
from constants import *


def solveroot(A, B, C):
    '''This function solves the root
    of the polynomial x^8+ax^6+bx^3+c=0 which 
    corresponds to the norm of r2, it needs to keep only
    the physical solutions (positive and not complex)'''
    coefs = np.array([1, 0, A, 0, 0, B, 0, 0, C])
    root = np.roots(coefs)
    ############################################
    ################TO DO######################
    ############################################


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

    # Creates the vector pointing to the station
    R1 = get_R(iodset[0].get_time(), lat1, h1)
    R2 = get_R(iodset[1].get_time(), lat2, h2)
    R3 = get_R(iodset[2].get_time(), lat3, h3)

    # Extracting observation time and time intervalls:
    t1 = iodset[0].get_time()
    t2 = iodset[1].get_time()
    t3 = iodset[2].get_time()
    tau = t3-t1
    tau1 = t1-t2
    tau3 = t3-t2

    # Unit vectors pointing to satellite
    e1 = iodset[0].get_e(lat1, lon1, R1)
    e2 = iodset[1].get_e(lat2, lon2, R2)
    e3 = iodset[2].get_e(lat3, lon3, R3)

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
    r2 = solveroot(a, b, c)

    rho1 = 1/D*((6*(D31*tau1/tau3+D21*tau/tau3)*r2**3+mu*D31 *
                (tau**2-tau1**2)*tau1/tau3)/(6*r2**3+mu*(tau**2-tau3**2))-D11)
    rho3 = 1/D*((6*(D13*tau3/tau1-D23*tau/tau1)*r2**3+mu*D13 *
                (tau**2-tau3**2)*tau3/tau1)/(6*r2**3+mu*(tau**2-tau3**2))-D33)

    # Vectors pointing to satellite
    r1 = rho1*e1+R1
    r2 = r2*e2+R2
    r3 = rho3*e3+R3

    return r1, r3


IODset = array(Create_IODs())
r1, r2 = find_r(IODset)
