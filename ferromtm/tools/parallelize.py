from multiprocessing import Pool
from functools import partial
from gridmap import grid_map


def parallel(function, partype="gridmap"):
    """
    Decorator to parallelize problems.

    Inputs
    ======
    function : the function that will be parallelized. The FIRST
        argument is the one to be iterated on (in parallel). The other
        arguments are the same in all the parallel runs of the function
        (they can be named or unnamedarguments).

    partype : the module used for parallelisation (gridmap, scoop or multiprocessing)

    Output
    ======
    A paralelized function. DO NOT NAME IT THE SAME AS THE INPUT
    FUNCTION.

    """
    if partype == "multiprocessing":
        # multiprocessing
        pool = Pool()
        parmap = pool.map
    else:
        # No parallelisation
        parmap = map

    def par(iterable_values, *args, **kwargs):
        mapfunc = partial(function, *args, **kwargs)
        result = [*parmap(mapfunc, iterable_values)]
        return result

    if partype == "gridmap":

        def par(iterable_values):
            result = grid_map(
                function,
                iterable_values,
                quiet=False,
                num_slots=1,
                temp_dir=u"./tmp",
                max_processes=len(iterable_values),
                queue="all.q",
                require_cluster=False,
                local=False,
                cleanup=True,
            )

            return result

    return par
