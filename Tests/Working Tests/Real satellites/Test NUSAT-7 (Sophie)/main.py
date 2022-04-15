import Gauss_algorithm
from IOD_class import *
from TLE_class import *
import math
import numpy as np
###############################################################


def main():
    # Creates three IOD
    IODset = np.array(Create_IODs())

    # Gauss orbit determination
    # precision
    prec = 0.00001
    # maximal iterations for Gauss algorithm
    itmax = 10000

    # Computes speed and position of satellite
    r2, v2 = Gauss_algorithm.find_r(IODset, itmax, prec)
    print("r2: ", math.sqrt(np.dot(r2, r2))-6361.54)
    print("v2: ", math.sqrt(np.dot(v2, v2)))

    # Corresponding observatioPositionn time
    t2 = IODset[1].get_time()
    TLEinput = [r2*1000, v2*1000, t2]

    # Create TLE
    tle = TLE(TLEinput)
    tle_list = tle.list_TLE()
    print(tle_list)

    # Create TLE file
    tle.TLE_format()


if __name__ == "__main__":
    main()
