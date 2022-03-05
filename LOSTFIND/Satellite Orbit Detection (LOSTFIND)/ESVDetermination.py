from IOD_class import *
from numpy import *


def find_r():
    '''This code uses three IOD' and then computes 
    two vectors r1 and r2 pointing to the satellite'''
    iodset = Create_IODs()

    # First estimations of n1 and n2:
    t1 = iodset(0).get_time_seconds()
    t2 = iodset(1).get_time_seconds()
    t3 = iodset(2).get_time_seconds()
    n1 = (t3-t2)/(t3-t1)
    n2 = (t2-t1)/(t3-t1)

    r1 = array([0, 0, 0])
    r2 = array([0, 0, 0])
    return r1, r2
