from IOD_class import *
from numpy import *
import os
import ephem
import datetime  
import math
import numpy as np


def get_R(epoch_time):
    '''This function reads the telescopes GPS position
    and then finds the coordinates in the absolute coordinate system'''
    script_dir = os.path.dirname(__file__)  # Location of python script
    rel_path = "LOCATION.txt"  # Name of .txt
    abs_file_path = os.path.join(script_dir, rel_path)
    Data = open(abs_file_path, "r")
    Station_Coordinates = Data.readlines()
    ################################################################
    ########################Rick part, may be wrong###################################
    ################################################################
    Re = 6378137  # Equatorial earth radius in meter
    f = 0.003353  # oblateness
    lat = Station_Coordinates[0]  # read latitudes GPS coordinates
    H = Station_Coordinates[2]  # read altitudes from GPS or google map 
    station =ephem.Observer()
    date_time = datetime.datetime.fromtimestamp( epoch_time )  
    station.date = date_time
    time = station.sidereal_time()
    R = ((Re/(math.sqrt(1-(2*f-f*f)*math.sin(lat)*math.sin(lat)))+H)*math.cos(lat)*math.cos(time),
         (Re/(math.sqrt(1-(2*f-f*f)*math.sin(lat)*math.sin(lat)))+H) *
         math.cos(lat)*math.sin(time),
         (Re*(1-f)*(1-f)/(math.sqrt(1-(2*f-f*f)*math.sin(lat)*math.sin(lat)))+H)*math.sin(lat))
    ################################################################
    #############################################################
    ###############################################################
    Data.close()
    return R


def find_rho(N1, N3, E1, E2, E3):
    '''This method takes the unit vectors and n1,n3 and computes
    the distances rho1,rho2,rho3 from the surface to the satellite'''
    RHO1 = 0
    RHO2 = 0
    RHO3 = 0
    return RHO1, RHO2, RHO3


def find_r(iodset):
    '''This code takes an array of IOD' and then computes 
    two vectors r1 and r2 pointing to the satellite'''

    # Creates the vector pointing to the station
    #R1 = get_R(iodset[0].get_time_unix())
    #R2 = get_R(iodset[1].get_time_unix())
    #R3 = get_R(iodset[2].get_time_unix())

    # First estimations of n1 and n2:
    #t1 = iodset[0].get_time_unix()
    #t2 = iodset[1].get_time_unix()
    #t3 = iodset[2].get_time_unix()

    ############################################################################
    #########HRemove before flight!!!!!!!!!##########
    ############################################################################

    t1 = 0
    t3 = 1
    t2 = 3
    R1 = array([1, 2, 3])
    R2 = array([1, 4, 3])
    R3 = array([1, 2, 5])
    # Unit vectors pointing to satellite
    e1 = iodset[0].get_e()
    e2 = iodset[1].get_e()
    e3 = iodset[2].get_e()

    # Usefull definitions
    d1 = cross(e2, e3)
    d2 = cross(e3, e1)
    d3 = cross(e1, e2)
    D = dot(e1, d1)

    if abs(t3-t1) and abs(t2-t1) > 0.00001:
        n1 = (t3-t2)/(t3-t1)
        n3 = (t2-t1)/(t3-t1)
    else:
        raise ValueError(
            "The IOD' must coorespond to different observation times")

    ############################################################################
    #########HERE THE CODE SHOULD DO A LOOP UNTIL PRECISION IS REACHED##########
    ############################################################################
    rho1, rho2, rho3 = find_rho(
        n1, n3, e1, e2, e3)
    # Computes the vector from earth center to satellite
    r1 = R1+rho1*e1
    r2 = R2+rho2*e2
    r3 = R3+rho3*e3
    # Next guess for n1, n3
    r1crossr3 = np.linalg.norm(cross(r1, r3), ord=2)
    r2crossr3 = np.linalg.norm(cross(r2, r3), ord=2)
    r1crossr2 = np.linalg.norm(cross(r1, r2), ord=2)
    n1 = r2crossr3/r1crossr3
    n3 = r1crossr2/r1crossr3

    return r1, r3


IODset = array(Create_IODs())
r1, r2 = find_r(IODset)
