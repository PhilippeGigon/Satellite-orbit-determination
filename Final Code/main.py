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
    prec = 0.0001
    # maximal iterations for Gauss algorithm
    itmax = 1000

    # Computes speed and position of satellite
    r2, v2 = Gauss_algorithm.find_r(IODset, itmax, prec, True)

    # Corresponding observatioPositionn time
    t2 = IODset[1].get_time()
    TLEinput = [r2, v2, t2]

    # Create TLE
    tle = TLE(TLEinput)

    # Create TLE file
    tle.TLE_format()


if __name__ == "__main__":
    main()
