from coupled2D import *
from ferromtm.tools.parallelize import *


if __name__ == "__main__":
    tstart = time.time()
    print("TSART: {}".format(tstart))
    out_cv = main_circle_conv(params[104])
    tend = time.time()
    et = tend - tstart
    print("TEND: {}".format(tend))
    print("ELAPSED: {}s".format(et))
