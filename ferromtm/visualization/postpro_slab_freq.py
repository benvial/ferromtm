from ferromtm.visualization.plots import *
from ferromtm.models.slab_freq import *

if __name__ == "__main__":

    arch = np.load(save_arch)
    out_meta = arch["out_meta"]
    out_meta_uncpl = arch["out_meta_uncpl"]

    R_hom = main_hom(np.array(fnorm), coupled=True)
    R_hom_uncpl = main_hom(np.array(fnorm), coupled=False)

    R_meta, T_meta, Q_meta, B_meta = extract_effs(out_meta)
    R_meta_uncpl, T_meta_uncpl, Q_meta_uncpl, B_meta_uncpl = extract_effs(
        out_meta_uncpl
    )

    plt.close("all")
    # plt.figure()
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5, 3))
    fig.set_rasterized(True)
    ax.set_rasterized(True)
    plt.plot(fnorm, R_meta, "-r", alpha=0.5, lw=4, label="MTM coupled")
    plt.plot(fnorm, R_hom, "--r", label="hom. coupled")
    plt.plot(fnorm, R_meta_uncpl, "-b", alpha=0.5, lw=4, label="MTM uncoupled")
    plt.plot(fnorm, R_hom_uncpl, "--b", label="hom. uncoupled")
    plt.legend(ncol=2)
    plt.xlabel("Normalized frequency $d/\lambda$")
    plt.ylabel("Reflection")
    plt.xlim((fnorm[0], 0.06))
    plt.ylim((0, 1.09))
    plt.tight_layout()
    fig.savefig("Rslab_freq.eps", rasterized=True, dpi=300)
