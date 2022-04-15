import Gauss_algorithm
import Laplace_algorithm
from IOD_class import *
from TLE_class import *
import math
import numpy as np


def main():
    prec = 0.0001  # precision
    itmax = 1000  # maximal iterations
    # Creates three IOD
    IODset = np.array(Create_IODs())
    print("Exercice 5.28, r=9759.8, v=6.0713, e=0.1, i=30")

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
    print('v2 Laplace:', math.sqrt(np.dot(v2, v2)))


if __name__ == "__main__":
    main()
