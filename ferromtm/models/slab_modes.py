from slab import *

maps_dir_ = "circ_rods"
maps_dir = os.path.join(data_folder, maps_dir_)


def load_arch(E_bias, coupling=True):
    f = 0.5
    fname = "maps_circle_f_{:.2f}_E_{:.2f}".format(f, E_bias)
    if not coupling:
        fname += "_uncoupled"
    filename = os.path.join(maps_dir, fname + ".npz")
    arch = np.load(filename)
    eps_hom = arch["eps_hom"]
    epsi_map = arch["epsi_map"]
    return epsi_map, eps_hom


fem.theta_deg = 0
fem.h_sup = 1  #: flt: "thickness" superstrate
fem.h_sub = 1  #: flt: "thickness" substrate


def main_modes(E_bias, coupled=True):
    epsi_map, eps_hom = load_arch(E_bias, coupled)
    fem.tmp_dir = tempfile.mkdtemp(prefix="/tmp/benjaminv.")
    fem.lambda0 = 150
    fem.lambda0search = 100
    fem.neig = 50

    ct = np.cos(fem.theta_deg * pi / 180)
    fem.h_pmltop = fem.lambda0 / ct  #: flt: thickness pml top
    fem.h_pmlbot = fem.lambda0 / ct  #: flt: thickness pml bot

    fem.analysis = "modal"
    fem.initialize()
    mesh = fem.make_mesh()
    make_epsi(fem, epsi_map, nincly)
    effs = []
    fem.compute_solution()
    ev = fem.postpro_eigenvalues()
    # fem.postpro_fields(filetype="pos")
    # fem.open_gmsh_gui()
    print("eigenvalues", ev)
    return ev


def main_meta(E_bias):
    return main_modes(E_bias)


def main_meta_uncpl(E_bias):
    return main_modes(E_bias, coupled=False)


# def main_hom(fnorm, coupled=True):
#     if coupled:
#         iter = 12  # periodic
#     else:
#         iter = 0
#     _, _, eps_hom = load_arch(iter)
#     lambda0 = fem.d/fnorm
#     r = rslab(fem.h_des, lambda0, fem.theta_deg, eps_hom[0,0], eps_hom[1,1])
#     return np.abs(r)**2


if __name__ == "__main__":

    run = False

    save_dir_ = "meta"
    save_arch = os.path.join(data_folder, save_dir_, "modes.npz")

    nE = 21
    Ebias = list(np.linspace(0, 2, nE))

    if run:
        from aotomat.tools.parallelize import *

        main_meta_par = parallel(main_meta, partype="gridmap")
        main_meta_uncpl_par = parallel(main_meta_uncpl, partype="gridmap")

        ev_meta = main_meta_par(Ebias)
        ev_meta_uncpl = main_meta_uncpl_par(Ebias)

        np.savez(save_arch, ev_meta=ev_meta, ev_meta_uncpl=ev_meta_uncpl)
    else:
        arch = np.load(save_arch)
        ev_meta = arch["ev_meta"]
        ev_meta_uncpl = arch["ev_meta_uncpl"]

    omega0 = 2 * pi / fem.d
    ev_meta_norma = np.array(ev_meta).ravel() / omega0
    ev_meta_uncpl_norma = np.array(ev_meta_uncpl).ravel() / omega0

    plt.figure()
    plt.plot(
        ev_meta_norma.real,
        ev_meta_norma.imag,
        "ro",
        alpha=0.5,
        lw=4,
        label="MTM coupled",
    )
    plt.plot(
        ev_meta_uncpl_norma.real,
        ev_meta_uncpl_norma.imag,
        "bo",
        alpha=0.5,
        lw=4,
        label="MTM uncoupled",
    )
    # plt.plot(angle, R_hom, "--r", label="hom. coupled")
    # plt.plot(angle, R_meta_uncpl, "-b", alpha=0.5, lw=4, label="MTM uncoupled")
    # plt.plot(angle, R_hom_uncpl, "--b", label="hom. uncoupled")
    plt.legend()
    plt.xlabel(r"Re $\omega$")
    plt.ylabel(r"Im $\omega$")
    # plt.xlim((0, 90))
    # plt.ylim((0, 1))
