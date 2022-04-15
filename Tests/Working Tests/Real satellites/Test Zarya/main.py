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
    prec = 0.1
    # maximal iterations for Gauss algorithm
    itmax = 1000

    # Computes speed and position of satellite
    r2, v2 = Gauss_algorithm.find_r(IODset, itmax, prec)
    print("r2: ", math.sqrt(np.dot(r2, r2))-6367.131)
    # Corresponding observatioPositionn time
    t2 = IODset[1].get_time()
    # In SI Units
    TLEinput = [r2*1000, v2*1000, t2]

    # Create TLE
    tle = TLE(TLEinput)
    tle_list = tle.list_TLE()
    print(tle_list)

    # Create TLE file
    tle.TLE_format()


if __name__ == "__main__":
    main()
