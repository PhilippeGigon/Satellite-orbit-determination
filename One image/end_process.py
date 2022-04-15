from detect_position import time_beginning, time_middle, time_ending
from utils.read_param import *
from altaz_to_IOD import *
import os

def creation_IOD_files(args, satellites):
    time, agvert, aghoriz, bouss, stnb, ststatut, texpos, focal, pixsize = read_param("parametres.txt")
    pscale = pixelscale(pixsize, focal)
    j = 0
    nameIOD = 'IOD_sat_numero'
    IODs = []
    for sat in satellites:
        j = j+1
        for i in range(0,len(sat.rhos)):
            tb = time_beginning(sat.times[i])
            tm = time_middle(sat.times[i], texpos)
            te = time_ending(sat.times[i], texpos)
            pbx = sat.beginningsx[i]
            pby = sat.beginningsy[i]
            pmx = sat.middlesx[i]
            pmy = sat.middlesy[i]
            pex = sat.endingsx[i]
            pey = sat.endingsy[i]
            [elb, azb] = positionAZEL(pbx, pby, agvert, aghoriz, bouss, pscale)
            [elm, azm] = positionAZEL(pmx, pmy, agvert, aghoriz, bouss, pscale)
            [ele, aze] = positionAZEL(pex, pey, agvert, aghoriz, bouss, pscale)

            elbstring = "%06d"%(elb)
            azbstring = "%07d"%(azb)
            elmstring = "%06d"%(elm)
            azmstring = "%07d"%(azm)
            elestring = "%06d"%(ele)
            azestring = "%07d"%(aze)

            sgb = ""  # signe elevation
            if elb > 0: sgb = "+"
            else: sgb = "-"
            
            sgm = "" 
            if elm > 0: sgm = "+"
            else: sgm = "-"

            sge = ""  # signe elevation
            if ele > 0: sge = "+"
            else: sge = "-"

            fnameIODb = './Processed/IOD/' + nameIOD[:-4] + str(j) + 'b.txt'
            fnameIODm = './Processed/IOD/' + nameIOD[:-4] + str(j) + 'm.txt'
            fnameIODe = './Processed/IOD/' + nameIOD[:-4] + str(j) + 'e.txt'
            angleforme = "04"
            tuncert = "18"
            puncert = "18"
            objnb = "12345 98 123ABC"
            optbehave = 'B'
            vismagsign = '+'
            vismag = 'MMm'
            maguncert = 'Mm'
            flashper = 'SSSsss'

            iodb = conv_altaz_to_IOD(elbstring, azbstring, sgb, angleforme, puncert, tb, tuncert, objnb,\
                optbehave, vismagsign, vismag, maguncert, flashper, stnb, ststatut, fnameIODb)
            IODs.append(iodb)

            iodm = conv_altaz_to_IOD(elmstring, azmstring, sgm, angleforme, puncert, tm, tuncert, objnb, \
                optbehave, vismagsign, vismag, maguncert, flashper, stnb, ststatut, fnameIODm)
            IODs.append(iodm)

            iode = conv_altaz_to_IOD(elestring, azestring, sge, angleforme, puncert, te, tuncert, objnb, \
                optbehave, vismagsign, vismag, maguncert, flashper, stnb, ststatut, fnameIODe)
            IODs.append(iode)
    return IODs


def Run_elfind_on_sats(pathelfind):
    '''
    os.system('ls')
    path='./Processed/IOD'
    path = path[2:len(path)]
    for fname in os.listdir(path):
        cmd = 'elfind ' + fname + '/'
        print(cmd)
        os.system(cmd)
    '''
    path = "./Processed/IOD"
    listnames = os.listdir(path)
    path3 = os.getcwd()
    os.chdir(pathelfind)
    for fname in listnames:
        fullname = path3 + '/Processed/IOD/' + fname
        print("Cannot run elfind.exe on Mac. Compile the cpp file in a Mac friendly way if needed.")
        #('elfind ' + fullname)