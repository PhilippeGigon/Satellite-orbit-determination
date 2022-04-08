from IOD_class import *
from numpy import *
import os
import datetime
import math
import numpy as np
import constants
from StationPositionVector import *
from Coordinates_conversion import *
#############################################################
# First version of orbit determination according to book:
# Satellite Orbits by Oliver Montbruck, this method never worked
#############################################################


def find_rho(N1, N3, D, D11, D12, D13, D21, D22, D23, D31, D32, D33):
    '''This method takes the vector products of the
    unit vectors and n1,n3 and computes the distances
    rho1,rho2,rho3 from the surface to the satellite'''
    RHO1 = (-1.0/(N1*D))*(N1*D11-D12+N3*D13)
    RHO2 = (1.0/D)*(N1*D21-D22+N3*D23)
    RHO3 = (-1.0/(N3*D))*(N1*D31-D32+N3*D33)
    return RHO1, RHO2, RHO3


def find_r(iodset):
    '''This code takes an array of IOD' and then computes
    two vectors r1 and r2 pointing to the satellite'''

    # Extracting location of station
    script_dir = os.path.dirname(__file__)  # Location of python script
    rel_path = "LOCATION_ELFIND2.txt"  # Name of .txt
    abs_file_path = os.path.join(script_dir, rel_path)
    Data = open(abs_file_path, "r")
    Station_Coordinates = Data.readlines()
    Data.close()
    lat1 = float(Station_Coordinates[0])*constants.degtorad
    lon1 = float(Station_Coordinates[1])*constants.degtorad
    h1 = float(Station_Coordinates[2])
    lat2 = float(Station_Coordinates[3])*constants.degtorad
    lon2 = float(Station_Coordinates[4])*constants.degtorad
    h2 = float(Station_Coordinates[5])
    lat3 = float(Station_Coordinates[6])*constants.degtorad
    lon3 = float(Station_Coordinates[7])*constants.degtorad
    h3 = float(Station_Coordinates[8])

    # Extracting observation time and time intervalls:
    t1 = iodset[0].get_time()
    t2 = iodset[1].get_time()
    t3 = iodset[2].get_time()
    time_s1 = get_lst(lon1, t1)
    time_s2 = get_lst(lon2, t2)
    time_s3 = get_lst(lon3, t3)
    t1 = 0
    t2 = 118.10
    t3 = 237.58

    tau = t3-t1
    tau1 = t1-t2
    tau3 = t3-t2
    # test
    time_s1 = 44.506*constants.degtorad
    time_s2 = 45.0*constants.degtorad
    time_s3 = 45.499*constants.degtorad

    # Creates the vector pointing to the station
    R1 = get_R(time_s1, lat1, h1)
    R2 = get_R(time_s2, lat2, h2)
    R3 = get_R(time_s3, lat3, h3)

    # Unit vectors pointing to satellite
    e1 = iodset[0].get_e(lat1, time_s1)
    e2 = iodset[1].get_e(lat2, time_s2)
    e3 = iodset[2].get_e(lat3, time_s3)
    # print(e1)
    # print(e2)
    # print(e3)

    # Usefull definitions
    d1 = cross(e2, e3)
    d2 = cross(e3, e1)
    d3 = cross(e1, e2)
    D11 = dot(d1, R1)
    D12 = dot(d1, R2)
    D13 = dot(d1, R3)
    D21 = dot(d2, R1)
    D22 = dot(d2, R2)
    D23 = dot(d2, R3)
    D31 = dot(d3, R1)
    D32 = dot(d3, R2)
    D33 = dot(d3, R3)
    D = dot(e3, d3)
    # Intial guess for n1, n3
    if abs(t3-t1) and abs(t2-t1) > 0.00000001:
        n1 = (t3-t2)/(t3-t1)
        n3 = (t2-t1)/(t3-t1)
        #print("Initial n1,n2:", n1, n3)
        n1old = 0
        n3old = 0

    else:
        raise ValueError(
            "The IOD' must coorespond to different observation times")

    ############################################################################
    #########HERE THE CODE SHOULD DO A LOOP UNTIL PRECISION IS REACHED##########
    ############################################################################
    # Difference between n and previous n smaller than epsilon-->stop
    epsilon = 0.000000001
    itmax = 100  # Maximum number of iterations
    iteration = 0
    # Does the loop as long as the n1,n3 change significat
    while True:
        iteration += 1
        rho1, rho2, rho3 = find_rho(
            n1, n3, D, D11, D12, D13, D21, D22, D23, D31, D32, D33)
        # Computes the vector from earth center to satellite
        r1 = R1+rho1*e1
        r3 = R3+rho3*e3
        r2 = R2+rho2*e2
        # print(r1)
        # print(R1)
        # print(R2)
        # print(R3)
        # Next guess for n1, n3
        r1crossr3 = np.sqrt(dot(cross(r1, r3), cross(r1, r3)))
        r2crossr3 = np.sqrt(dot(cross(r2, r3), cross(r2, r3)))
        r1crossr2 = np.sqrt(dot(cross(r1, r2), cross(r1, r2)))
        print("===")
       # print(r1)
        # print(r3)
        print(r1crossr3)
        print(r2crossr3)
        print(float(r2crossr3/r1crossr3))
        print(float(r1crossr2/r1crossr3))
        print(n1)
        print(n3)
        print("===")
        if r1crossr3 == 0:
            raise ValueError(
                "Two vectors are parallel!")
        n1old = n1
        n3old = n3
        n1 = float(r2crossr3/r1crossr3)
        n3 = float(r1crossr2/r1crossr3)

        if abs(n1-n1old) < epsilon and abs(n3-n3old) < epsilon or iteration > itmax:
            print("Number of iterations: ", iteration,
                  "difference:", abs(n1-n1old))
            break

    r1 = R1+rho1*e1
    r3 = R3+rho3*e3
    return r1, r3
