import Gauss_algorithm
from IOD_class import *
from TLE_class import *
import math
import numpy as np
import Laplace_algorithm


def main():
    prec = 0.0001  # precision
    itmax = 1000  # maximal iterations
    # Creates three IOD
    IODset = np.array(Create_IODs())
    print("Exercice 5.22, r=6743, v=7.6922, i=52, e=0.001")

    # Computes two vectors pointing to the satellite with gauss method
    r2, v2 = Gauss_algorithm.find_r(
        IODset, itmax, prec)
    print('r2 Gauss:', math.sqrt(np.dot(r2, r2)))
    print('v2 Gauss:', math.sqrt(np.dot(v2, v2)))

    t2 = IODset[1].get_time()
    TLEinput = [r2, v2, t2]
    # Create TLE
    tle = TLE(TLEinput)
    # Create TLE file
    tle.TLE_format()

    r2, v2 = Laplace_algorithm.rv(IODset, 2)
    print('r2 Laplace:', math.sqrt(np.dot(r2, r2)))


if __name__ == "__main__":
    main()
