import datetime
import math
import numpy as np
import ephem

def get_R(epoch_time, lat, H):
    ################################################################
    ########################Rick part, may be wrong###################################
    ################################################################
    Re = 6378137  # Equatorial earth radius in meter
    f = 0.003353  # oblateness
    # read latitudes GPS coordinates

    station = ephem.Observer()
    date_time = datetime.datetime.fromtimestamp(epoch_time)
    station.date = date_time
    time = station.sidereal_time()
    R = ((Re/(math.sqrt(1-(2*f-f*f)*math.sin(lat)*math.sin(lat)))+H)*math.cos(lat)*math.cos(time),
         (Re/(math.sqrt(1-(2*f-f*f)*math.sin(lat)*math.sin(lat)))+H) *
         math.cos(lat)*math.sin(time),
         (Re*(1-f)*(1-f)/(math.sqrt(1-(2*f-f*f)*math.sin(lat)*math.sin(lat)))+H)*math.sin(lat))
    ################################################################
    #############################################################
    ###############################################################
    return R
