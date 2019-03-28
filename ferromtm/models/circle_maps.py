from ferromtm.models.coupled2D import *
from ferromtm.tools.parallelize import *


def main_circle(params, save=False, coupling=True):
    E_bias, f = params
    print("Parameters: E = {:.2f}MV/m - f = {:.2f} ".format(E_bias, f))
    eps_hom, epsi, E, fem_hom, fem_es, epsi_map = main(
        f, E_bias, coupling=coupling, postmaps=True
    )

    if save:
        fname = "maps_circle_f_{:.2f}_E_{:.2f}".format(f, E_bias)
        if not coupling:
            fname += "_uncoupled"
        saveddir = os.path.join(data_folder, "circ_rods")
        try:
            os.mkdir(saveddir)
        except FileExistsError:
            pass
        np.savez(
            os.path.join(saveddir, fname + ".npz"),
            eps_hom=eps_hom,
            f=f,
            E_bias=E_bias,
            epsi=epsi,
            epsi_map=epsi_map,
            E=E,
        )
    return eps_hom, epsi, E, fem_hom, fem_es


def main_per_uncoupled(params):
    return main_circle(params, coupling=False, save=True)


def main_per_coupled(params):
    return main_circle(params, save=True, coupling=True)


nE = 21
Ebias = np.linspace(0, 2, nE)
nF = 1
Fincl = 0.5
E1, F1 = np.meshgrid(Ebias, Fincl)
params = np.vstack((E1.ravel(), F1.ravel())).T


if __name__ == "__main__":

    run_parallel_coupled = parallel(main_per_coupled, partype="gridmap")
    run_parallel_uncoupled = parallel(main_per_uncoupled, partype="gridmap")
    tstart = time.time()
    print("TSART: {}".format(tstart))
    out_coupled = run_parallel_coupled(params)
    out_uncoupled = run_parallel_uncoupled(params)
    tend = time.time()
    et = tend - tstart
    print("TEND: {}".format(tend))
    print("ELAPSED: {}s".format(et))
