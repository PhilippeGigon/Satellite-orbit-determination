from ESVDetermination import *
from IOD_class import *


def main():
    # Creates three IOD
    IODset = array(Create_IODs())
    # Computes two vectors pointing to the satellite
    r1, r2 = find_r(IODset)
    # Corresponding observation time
    t1 = IODset[0].get_time()
    t2 = IODset[2].get_time()
    Position = [r1, r2, t1, t2]
    print(sqrt(dot(r2, r2))/1000.0-6371.0)


if __name__ == "__main__":
    main()
