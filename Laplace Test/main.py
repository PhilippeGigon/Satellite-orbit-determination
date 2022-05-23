import Gauss_algorithm
import Laplace_algorithm
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
    r2G, v2G = Gauss_algorithm.find_r(IODset, itmax, prec, True)
    print("Gauss r: ", math.sqrt(np.dot(r2G, r2G))/1000)
    print("Gauss v: ", math.sqrt(np.dot(v2G, v2G))/1000)

    r2L, v2L = Laplace_algorithm.rv(IODset, 2)
    print("Laplace r: ", math.sqrt(np.dot(r2L, r2L)))
    print("Laplace v: ", math.sqrt(np.dot(v2L, v2L)))

    print("rgauss-rlaplace: ", math.sqrt(np.dot(r2L-r2G/1000, r2L-r2G/1000)))


if __name__ == "__main__":
    main()
