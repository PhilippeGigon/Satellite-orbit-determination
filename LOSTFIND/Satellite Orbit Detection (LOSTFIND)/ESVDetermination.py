from IOD_class import *
from numpy import *
import os


def get_R():
    '''This function reads the telescopes GPS position
    and then finds the coordinates in the absolute coordinate system'''
    ################################################################
    ########################TO DO###################################
    ################################################################
    script_dir = os.path.dirname(__file__)  # Location of python script
    rel_path = "TelescopeLocation.txt"  # Name of .txt
    abs_file_path = os.path.join(script_dir, rel_path)
    Data = open(abs_file_path, "r")
    Station_Coordinates = Data.readlines()
    Data.close()


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
