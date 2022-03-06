from ESVDetermination import *
from IOD_class import *


def main():
    # Creates three IOD
    IODset = Create_IODs()

    # Computes two vectors pointing to the satellite
    r1, r2 = find_r(IODset)


if __name__ == "__main__":
    main()
