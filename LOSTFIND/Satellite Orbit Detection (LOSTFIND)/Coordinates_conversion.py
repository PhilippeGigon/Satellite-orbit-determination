import math
import ephem

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
    az = math.atan(-math.cos(dec)*math.cos(lat)*math.sin(H)/(math.sin(dec)-math.sin(lat)*math.sin(alt)))  
    return az
    	
