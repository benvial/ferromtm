from slab import *


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


if __name__ == "__main__":

    run = False

    fnorm = list(np.linspace(1 / 150, 1 / 10, 150))
    save_dir_ = "meta"
    save_arch = os.path.join(data_folder, save_dir_, "efficiencies_freq.npz")

    if run:
        from aotomat.tools.parallelize import *

        main_meta_par = parallel(main_meta, partype="gridmap")
        main_meta_uncpl_par = parallel(main_meta_uncpl, partype="gridmap")

        out_meta = main_meta_par(fnorm)
        out_meta_uncpl = main_meta_uncpl_par(fnorm)

        np.savez(save_arch, out_meta=out_meta, out_meta_uncpl=out_meta_uncpl)
    else:
        arch = np.load(save_arch)
        out_meta = arch["out_meta"]
        out_meta_uncpl = arch["out_meta_uncpl"]

    R_hom = main_hom(np.array(fnorm), coupled=True)
    R_hom_uncpl = main_hom(np.array(fnorm), coupled=False)

    R_meta, T_meta, Q_meta, B_meta = extract_effs(out_meta)
    R_meta_uncpl, T_meta_uncpl, Q_meta_uncpl, B_meta_uncpl = extract_effs(
        out_meta_uncpl
    )

    plt.figure()
    plt.plot(fnorm, R_meta, "-r", alpha=0.5, lw=4, label="MTM coupled")
    plt.plot(fnorm, R_hom, "--r", label="hom. coupled")
    plt.plot(fnorm, R_meta_uncpl, "-b", alpha=0.5, lw=4, label="MTM uncoupled")
    plt.plot(fnorm, R_hom_uncpl, "--b", label="hom. uncoupled")
    plt.legend()
    plt.xlabel("Normalized frequency $d/\lambda$")
    plt.ylabel("Reflection")
    # plt.xlim((0, 90))
    plt.ylim((0, 1))
