from ferromtm.visualization.plots import *
from ferromtm.models.slab_modes import *

if __name__ == "__main__":
    arch = np.load(save_arch)
    ev_meta = arch["ev_meta"]
    ev_meta_uncpl = arch["ev_meta_uncpl"]

    omega0 = 2 * pi / fem.d
    ev_meta_norma = np.array(ev_meta).conj() / omega0
    ev_meta_uncpl_norma = np.array(ev_meta_uncpl).conj() / omega0

    ev_slab_uncpl = [main_modes_analytical(e, coupled=False) for e in Ebias]
    ev_slab_uncpl_norma = np.array(ev_slab_uncpl) / omega0
    ev_slab = [main_modes_analytical(e, coupled=True) for e in Ebias]
    ev_slab_norma = np.array(ev_slab) / omega0

    dx = +0.002
    x1, y1 = ev_slab_norma.real[13][1] + dx, ev_slab_norma.imag[13][1]
    x2, y2 = ev_slab_norma.real[5][1] + dx, ev_slab_norma.imag[5][1]

    locs = x1, y1, x2, y2

    def draw_arrow(ax, locs, connectionstyle, label=None):
        x1, y1, x2, y2 = locs
        ax.annotate(
            "",
            xy=(x1, y1),
            xycoords="data",
            xytext=(x2, y2),
            textcoords="data",
            arrowprops=dict(
                arrowstyle="->",
                color="0.0",
                shrinkA=5,
                shrinkB=5,
                patchA=None,
                patchB=None,
                connectionstyle=connectionstyle,
            ),
        )
        ax.annotate(
            label,
            xy=(x1, y1),
            xycoords="data",
            xytext=((x1 + x2) / 2, (y1 + y2) / 2),
            textcoords="data",
        )

    plt.close("all")
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(5, 3))
    plt.plot(
        ev_meta_norma.real.ravel(),
        ev_meta_norma.imag.ravel(),
        "ko",
        alpha=0.3,
        lw=4,
        label="MTM coupled",
    )
    plt.plot(
        ev_slab_norma.real.ravel(),
        ev_slab_norma.imag.ravel(),
        "rx",
        alpha=1,
        lw=4,
        label="slab coupled",
    )
    plt.xlim((0.01, 0.06))
    plt.ylim((-0.004, 0.0))
    draw_arrow(ax, locs, "arc3,rad=0.", label=r"$E_0 \nearrow$")
    # plt.plot(
    #     ev_meta_uncpl_norma.real.ravel(),
    #     ev_meta_uncpl_norma.imag.ravel(),
    #     "bo",
    #     alpha=0.5,
    #     lw=4,
    #     label="MTM uncoupled",
    # )
    # plt.plot(
    #     ev_slab_uncpl_norma.real.ravel(),
    #     ev_slab_uncpl_norma.imag.ravel(),
    #     "bx",
    #     alpha=1,
    #     lw=4,
    #     label="slab uncoupled",
    # )
    # plt.plot(angle, R_hom, "--r", label="hom. coupled")
    # plt.plot(angle, R_meta_uncpl, "-b", alpha=0.5, lw=4, label="MTM uncoupled")
    # plt.plot(angle, R_hom_uncpl, "--b", label="hom. uncoupled")
    plt.legend()
    plt.xlabel(r"Re $\omega/\omega_0$")
    plt.ylabel(r"Im $\omega/\omega_0$")
    # plt.xlim((0, 90))
    # plt.ylim((0, 1))

    fig.set_rasterized(True)
    ax.set_rasterized(True)
    plt.tight_layout()
    fig.savefig("qnms_slab.eps", rasterized=True, dpi=300)
