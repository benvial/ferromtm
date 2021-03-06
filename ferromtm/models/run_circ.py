from ferromtm.models.coupled2D import *
from ferromtm.tools.parallelize import *


def main_per_uncoupled(params):
    return main_circle(params, coupling=False, save=True)


def main_per_coupled(params):
    return main_circle(params, save=True, coupling=True)


def test(params):
    return 0


if __name__ == "__main__":

    run_parallel_coupled = parallel(main_per_coupled, partype="gridmap")
    run_parallel_uncoupled = parallel(main_per_uncoupled, partype="gridmap")
    test_parallel = parallel(test)
    tstart = time.time()
    print("TSART: {}".format(tstart))
    # out_coupled = run_parallel_coupled(params)
    # out_uncoupled = run_parallel_uncoupled(params)
    out_test = test_parallel(params)
    tend = time.time()
    et = tend - tstart
    print("TEND: {}".format(tend))
    print("ELAPSED: {}s".format(et))
