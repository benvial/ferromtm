from ferromtm.models.slab import *

fem.theta_deg = 20


def main_frequency(fnorm, coupled=True):
    # fnorm = period/lambda0
    if coupled:
        iter = 12  # periodic
    else:
        iter = 0
    epsi_map, E_map, eps_hom = load_arch(iter)
    fem.tmp_dir = tempfile.mkdtemp(prefix="/tmp/benjaminv.")
    fem.lambda0 = fem.d / fnorm

    ct = np.cos(fem.theta_deg * pi / 180)
    fem.h_pmltop = fem.lambda0 / ct  #: flt: thickness pml top
    fem.h_pmlbot = fem.lambda0 / ct  #: flt: thickness pml bot
    fem.initialize()
    mesh = fem.make_mesh()
    make_epsi(fem, epsi_map, nincly)
    effs = []
    fem.compute_solution()
    effs = fem.diffraction_efficiencies()
    print("efficiencies", effs)
    return effs


def main_meta(fnorm):
    return main_frequency(fnorm)


def main_meta_uncpl(fnorm):
    return main_frequency(fnorm, coupled=False)


def main_hom(fnorm, coupled=True):
    if coupled:
        iter = 12  # periodic
    else:
        iter = 0
    _, _, eps_hom = load_arch(iter)
    lambda0 = fem.d / fnorm
    r = rslab(fem.h_des, lambda0, fem.theta_deg, eps_hom[0, 0], eps_hom[1, 1])
    return np.abs(r) ** 2


fnorm = list(np.linspace(1 / 150, 1 / 10, 150))
save_dir_ = "meta"
save_arch = os.path.join(data_folder, save_dir_, "efficiencies_freq.npz")


if __name__ == "__main__":

    from ferromtm.tools.parallelize import *

    main_meta_par = parallel(main_meta, partype="gridmap")
    main_meta_uncpl_par = parallel(main_meta_uncpl, partype="gridmap")

    out_meta = main_meta_par(fnorm)
    out_meta_uncpl = main_meta_uncpl_par(fnorm)

    np.savez(save_arch, out_meta=out_meta, out_meta_uncpl=out_meta_uncpl)
