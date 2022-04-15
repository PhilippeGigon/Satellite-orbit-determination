import math
import ephem
import datetime

def get_lst(lon,epoch_time):
    station = ephem.Observer()
    date_time = datetime.datetime.fromtimestamp(epoch_time)
    station.date=date_time
    time=station.sidereal_time()+lon
    return time

def get_alt(dec,rad,lat,time): 
    H = time-rad
    alt = math.asin(math.sin(dec)*math.sin(lat)+math.cos(lat)*math.cos(dec)*math.cos(H))   
    return alt

def get_az(alt,dec,rad,lat,time): 
    H = time-rad
    az = math.acos(((math.sin(dec)-math.sin(lat)*math.sin(alt))/(math.cos(lat)*math.cos(alt)))) 
    if (math.sin(H)>0):
        az = 2*math.pi-az 
    return az    
    	         	
def get_dec(alt,az,lat): 
    dec = math.asin(math.sin(alt)*math.sin(lat)+math.cos(alt)*math.cos(lat)*math.cos(az))
    return dec

def get_rad(alt,az,lat,dec): 
    H = math.acos((math.sin(alt)-math.sin(lat)+math.cos(dec))/(math.cos(lat)*math.cos(dec)))
    rad = time-H
    return rad
    
    
