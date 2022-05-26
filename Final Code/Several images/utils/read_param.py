from detect_position import time_from_str_to_datetime

def read_param(param_filename):

    parametres = open(param_filename, "r")

    entete = parametres.readline().rstrip('\n\r').split(",")
    param1 = parametres.readline().rstrip('\n\r').split(",")
    param2 = parametres.readline().rstrip('\n\r').split(",")
    param3 = parametres.readline().rstrip('\n\r').split(",")

    boussoleidx = entete.index("Boussole")
    anglevertidx = entete.index("AngleVert")
    anglehorizidx = entete.index("AngleHoriz")
    focalidx = entete.index("Focale")
    pixelsizeidx = entete.index("Pixsize")
    stationnbidx = entete.index("StationNumber")
    stationstatutidx = entete.index("StationStatut")
    timeexposidx = entete.index("TmpsExposition") #in seconds
    timeidx = entete.index("Time")
    #nameinputidx=entete.index("Input")
    #nameoutputidx=entete.index("Output")

    focal = float(param1[focalidx])
    pixsize = float(param1[pixelsizeidx])
    agvert = float(param1[anglevertidx])
    aghoriz = float(param1[anglehorizidx])
    bouss = float(param1[boussoleidx])
    stnb = param1[stationnbidx]
    ststatut = param1[stationstatutidx]
    texpos = float(param1[timeexposidx])
    timestr1 = param1[timeidx]
    timestr2 = param2[timeidx]
    timestr3 = param3[timeidx]

    times = [timestr1, timestr2, timestr3]
    #inputname=param[nameinputidx]
    #outputname=param[nameoutputidx]

    return times,agvert,aghoriz,bouss,stnb,ststatut,texpos,focal,pixsize