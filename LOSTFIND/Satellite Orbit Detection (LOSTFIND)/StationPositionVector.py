import datetime
import math
import numpy as np
import ephem
import constants

def get_R(time, lat, H):
    ################################################################
    ########################Rick part, may be wrong###################################
    ################################################################
    # read latitudes GPS coordinates
    #local_time 
    R = ((constants.Re/(math.sqrt(1-(2*constants.f-constants.f*constants.f)*math.sin(lat)*math.sin(lat)))+H)*math.cos(lat)*math.cos(time),
    constants.Re/(math.sqrt(1-((2*constants.f-constants.f*constants.f)*math.sin(lat)*math.sin(lat)))+H) *
         math.cos(lat)*math.sin(time),
         (constants.Re*(1-constants.f)*(1-constants.f)/(math.sqrt(1-(2*constants.f-constants.f*constants.f)*math.sin(lat)*math.sin(lat)))+H)*math.sin(lat))
    ################################################################
    #############################################################
    ###############################################################
    return R
