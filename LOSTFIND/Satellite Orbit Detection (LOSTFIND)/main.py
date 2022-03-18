from ESVDetermination import *
from IOD_class import *
from TLE_class import *


def main():
    # Creates three IOD
    IODset = array(Create_IODs())
    # Computes two vectors pointing to the satellite
    r1, r2 = find_r(IODset)
    # Corresponding observation time
    t1 = IODset[0].get_time()
    t2 = IODset[2].get_time()
    Position = [r1, r2, t1, t2]
    print(sqrt(dot(r1, r1))/1000.0-6371.0)

    # Create TLE
    prec = 0.1  # precision
    tle = TLE(Position, prec)
    tle_list = tle.list_TLE()
    # print(tle_list)


if __name__ == "__main__":
    main()
