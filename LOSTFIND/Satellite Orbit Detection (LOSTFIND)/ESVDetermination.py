from IOD_class import *
from numpy import *
import os
import ephem
import math


def get_R(date_string):
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
    station = ephem.Observer()
    station.lon = Station_Coordinates[1]  # should be a str i think
    station.lat = Station_Coordinates[0]
    station.date = time
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


def find_r(iodset):
    '''This code takes an array of IOD' and then computes 
    two vectors r1 and r2 pointing to the satellite'''

    # First estimations of n1 and n2:
    t1 = iodset[0].get_time_seconds()
    t2 = iodset[1].get_time_seconds()
    t3 = iodset[2].get_time_seconds()
    print(t1, t2, t3)
    if abs(t3-t1) and abs(t2-t1) > 0.00001:
        n1 = (t3-t2)/(t3-t1)
        n2 = (t2-t1)/(t3-t1)
    else:
        raise ValueError(
            "The IOD' must coorespond to different observation times")

    r1 = array([0, 0, 0])
    r2 = array([0, 0, 0])
    return r1, r2
