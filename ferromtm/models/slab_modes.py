from ferromtm.models.slab import *


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


def main_modes_analytical(E_bias, coupled=True):
    _, eps_hom = load_arch(E_bias, coupled)
    h, epsxx, epsyy = fem.h_des, np.conj(eps_hom[0, 0]), np.conj(eps_hom[1, 1])
    n = 7
    n_ = np.array(range(7))
    alpha = (np.sqrt(epsxx) + 1) / (np.sqrt(epsxx) - 1)
    ev = (n_ * pi - 1j * np.log(alpha)) / (h * np.sqrt(epsxx))
    return ev


save_dir_ = "meta"
save_arch = os.path.join(data_folder, save_dir_, "modes.npz")

nE = 21
Ebias = list(np.linspace(0, 2, nE))

if __name__ == "__main__":
    from ferromtm.tools.parallelize import *

    main_meta_par = parallel(main_meta, partype="gridmap")
    main_meta_uncpl_par = parallel(main_meta_uncpl, partype="gridmap")

    ev_meta = main_meta_par(Ebias)
    ev_meta_uncpl = main_meta_uncpl_par(Ebias)

    np.savez(save_arch, ev_meta=ev_meta, ev_meta_uncpl=ev_meta_uncpl)
