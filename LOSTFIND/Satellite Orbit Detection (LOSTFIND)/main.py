import SatelitePositionDetermination_Gauss
from IOD_class import *
from TLE_class import *
import math
import numpy as np


def main():
    # Creates three IOD
    IODset = np.array(Create_IODs())

    # Computes two vectors pointing to the satellite with gauss method
    r1, r2 = SatelitePositionDetermination_Gauss.find_r(IODset)

    # Corresponding observation time
    t1 = IODset[0].get_time()
    t2 = IODset[2].get_time()

    #Position = [r1, r2, t1, t2]
    # Create TLE
    # prec = 0.1  # precision
    #tle = TLE(Position, prec)
    #tle_list = tle.list_TLE()
    # print(tle_list)

    # Create TLE file
    # tle.TLE_format()


if __name__ == "__main__":
    main()
