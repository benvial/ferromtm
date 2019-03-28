from ferromtm.models.coupled2D import *
from ferromtm.tools.parallelize import *
import sys

if __name__ == "__main__":
    if len(sys.argv) > 1:
        percase = sys.argv[1]
    else:
        percase = "per"

    tstart = time.time()
    if percase == "rand":
        out_cv = main_random_conv(params[104])
    else:
        out_cv = main_circle_conv(params[104])
    tend = time.time()
    et = tend - tstart
    print("ELAPSED: {}s".format(et))
