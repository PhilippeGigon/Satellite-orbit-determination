import os
from datetime import datetime

#   Contains functions to create param_img from the name of the folder containing all the images

def create_param_img(times, foldername): #time=real time ?
    param_img = open("param_img.txt","w+")
    foldername = '../' + foldername  #
    filenames = sorted([f for f in os.listdir(foldername) if f.endswith('.jpg')])
    for j in range(0, len(filenames)):
        time_string = str(times[j])
        param_img.write(filenames[j] + "    " + time_string + '\n')
    param_img.close()

def get_rasp_time(filename):
    print('filename:', filename)
    time = 0
    parameters = open('parameters.txt', "r" )
    return time
