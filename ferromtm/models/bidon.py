#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT
# import logging
# from gridmap import grid_map


# logging.captureWarnings(True)
# logging.basicConfig(
#     format=("%(asctime)s - %(name)s - %(levelname)s - " + "%(message)s"),
#     level=logging.DEBUG,
# )

import numpy as np
import time


def main_rand(params, save=False, coupling=True):
    print("fuck")
    return params


def main_rand_coupled(params):
    test = main_rand(params, save=True, coupling=True)
    return test


# function_parallel = main_rand_coupled

nE = 11
Ebias = np.linspace(0, 2, nE)
nF = 5
Fincl = np.linspace(0.1, 0.5, nF)
E1, F1 = np.meshgrid(Ebias, Fincl)
params = np.vstack((E1.ravel(), F1.ravel())).T
from aotomat.tools.parallelize import *

# def run_parallel(params):
#     result = grid_map(
#         main_rand_coupled,
#         params,
#         quiet=False,
#         num_slots=1,
#         temp_dir=u"./tmp",
#         max_processes=len(params),
#         queue="all.q",
#         require_cluster=True,
#         local=True,
#         cleanup=True,
#     )
#     return result


if __name__ == "__main__":
    run_parallel = parallel(main_rand_coupled, partype="gridmap")

    # f_parallel(params)
    tstart = time.time()
    print("TSART: {}".format(tstart))
    out = run_parallel(params)
    print(out)
    # f_parallel(params, save=True, coupling=False)
    tend = time.time()
    et = tend - tstart
    print("TEND: {}".format(tend))
    print("ELAPSED: {}s".format(et))
