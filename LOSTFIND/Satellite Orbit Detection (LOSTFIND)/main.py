from ESVDetermination import *
from IOD_class import *


def main():
    # Creates three IOD
    IODset = Create_IODs()
    for k in range(0, 3):
        R(k) = get_R(IODset(k).get_stringdate())
    # Computes two vectors pointing to the satellite
    #r1, r2 = find_r(IODset)


if __name__ == "__main__":
    main()
