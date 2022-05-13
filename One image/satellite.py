import numpy as np
import datetime
from detect_position import time_from_str_to_datetime, time_from_datetime_to_str
import math

def norm(trace_a, trace_b):
    (rho, theta, time) = trace_a
    (rhobis, thetabis, timebis) = trace_b
    #   provide rho_char, theta_char and time_char
    rho_char = 1
    theta_char = 1e-2
    time = time_from_str_to_datetime(time)
    timebis = time_from_str_to_datetime(timebis)
    time_char = datetime.timedelta(minutes = 4)
    #   end
    delta = abs(time - timebis)
    if delta > time_char:
        return 1e6
    else:
        print('NORM',math.sqrt((rho-rhobis)/rho_char*(rho-rhobis)/rho_char+(theta-thetabis)/theta_char*(theta-thetabis)/theta_char))
        return math.sqrt((rho-rhobis)/rho_char*(rho-rhobis)/rho_char+(theta-thetabis)/theta_char*(theta-thetabis)/theta_char)

class Satellite:
    def __init__(self, trace, beginning, middle, ending):
        (rhoss, thetass, timess) = trace
        (beginningsxx, beginningsyy) = beginning
        (middlesxx,middlesyy) = middle
        (endingsxx, endingsyy) = ending
        self.rhos = [rhoss]
        self.thetas = [thetass]
        self.times = [timess]
        
        self.beginningsx = []
        self.beginningsx.append(beginningsxx)
        self.beginningsy = []
        self.beginningsy.append(beginningsyy)

        self.middlesx = []
        self.middlesx.append(middlesxx)
        self.middlesy = []
        self.middlesy.append(middlesyy)

        self.endingsx = []
        self.endingsx.append(endingsxx)
        self.endingsy = []
        self.endingsy.append(endingsyy)


    def same_satellite(self, trace, beginning, middle, ending, threshold = 1):
        #   can modify to take average over the satellite's traces
        (Beginningsxx, Beginningsyy) = beginning
        (Middlesxx, Middlesyy) = middle
        (Endingxx, Endingyy) = ending
        rhos = self.rhos
        rho = rhos[0]
        thetas = self.thetas
        theta = thetas[0]
        times = self.times
        time = times[0]
        same = (norm(trace, (rho, theta,time)) < threshold)
        if same:
            self.rhos.append(rho)
            self.thetas.append(theta)
            self.times.append(time)
            self.beginningsx.append(Beginningsxx)
            self.beginningsy.append(Beginningsyy)
            self.middlesx.append(Middlesxx)
            self.middlesy.append(Middlesyy)
            self.endingsx.append(Endingxx)
            self.endingsy.append(Endingyy)
        return same