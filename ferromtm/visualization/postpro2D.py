from ferromtm.models.coupled2D import *
from ferromtm.visualization.plots import *


n_bulk = tunability(epsilonr_ferroelectric(Ebias), epsilonr_ferroelectric(0))
K_bulk = cqf(epsilonr_ferroelectric(Ebias), epsilonr_ferroelectric(0))
eps_00 = epsilonr_ferroelectric(0).real
tandelta_00 = -epsilonr_ferroelectric(0).imag / epsilonr_ferroelectric(0).real


def load_arch(E_bias, f, coupling=True):
    name = "circle_f_{:.2f}_E_{:.2f}".format(f, E_bias)
    if not coupling:
        name += "_uncoupled"
    saveddir = os.path.join(data_folder, "circ_rods")
    # saveddir = os.path.join(data_folder, "circ_rods_pattern")
    filename = os.path.join(saveddir, name + ".npz")
    arch = np.load(filename)
    eps_hom = arch["eps_hom"]
    epsi = arch["epsi"]
    E = arch["E"]
    return epsi, E, eps_hom


iplot = [0, 1, 2, 3, 4]


def load_results():
    Eps_eff_lin = np.zeros((nE, nF), dtype=complex)
    Eps_eff_nl = np.zeros((nE, nF), dtype=complex)
    Aniso_factor_lin = np.zeros((nE, nF), dtype=float)
    Aniso_factor_nl = np.zeros((nE, nF), dtype=float)
    #

    for iR, f in enumerate(Fincl):
        for iE, E_bias in enumerate(Ebias):
            _, _, eps_hom = load_arch(E_bias, f, coupling=True)
            eps_eff_nl = eps_hom[0, 0]
            Eps_eff_nl[iE, iR] = eps_eff_nl
            Aniso_factor_nl[iE, iR] = aniso_factor(eps_hom)

            _, _, eps_hom_lin = load_arch(E_bias, f, coupling=False)
            eps_eff_lin = eps_hom_lin[0, 0]
            Eps_eff_lin[iE, iR] = eps_eff_lin
            Aniso_factor_lin[iE, iR] = aniso_factor(eps_hom_lin)

    Eps_eff_lin = Eps_eff_lin[:, iplot]
    Eps_eff_nl = Eps_eff_nl[:, iplot]
    Aniso_factor_lin = Aniso_factor_lin[:, iplot]
    Aniso_factor_nl = Aniso_factor_nl[:, iplot]

    norm_eps = Eps_eff_lin.real / eps_00
    norm_loss = -Eps_eff_lin.imag / Eps_eff_lin.real / tandelta_00
    norm_eps_nl = Eps_eff_nl.real / eps_00
    norm_loss_nl = -Eps_eff_nl.imag / Eps_eff_nl.real / tandelta_00
    n_lin = tunability(Eps_eff_lin, Eps_eff_lin[0, :])
    n_nl = tunability(Eps_eff_nl, Eps_eff_nl[0, :])
    n_norm_lin = np.array([n / n_bulk for n in n_lin.T]).T
    n_norm_nl = np.array([n / n_bulk for n in n_nl.T]).T

    K_lin = cqf(Eps_eff_lin, Eps_eff_lin[0, :])
    K_nl = cqf(Eps_eff_nl, Eps_eff_nl[0, :])

    Knorm_lin = np.array([(k) for k in K_lin.T]).T
    Knorm_nl = np.array([(k) for k in K_nl.T]).T
    # K_nl[K_nl==0]=1e-12
    return (
        norm_eps,
        norm_eps_nl,
        norm_loss,
        norm_loss_nl,
        n_norm_lin,
        n_norm_nl,
        Aniso_factor_lin,
        Aniso_factor_nl,
    )


def plot_eff_par(fig, ax):
    ax[0][0].plot(Ebias, norm_eps, ls="--")
    ax[0][0].set_prop_cycle(None)
    ax[0][0].plot(Ebias, norm_eps_nl, ls="-")
    ax[1][0].plot(Ebias, norm_loss, ls="--")
    ax[1][0].set_prop_cycle(None)
    ax[1][0].plot(Ebias, norm_loss_nl, ls="-")

    ax[1][0].set_xlabel("bias electric field (MV/m)")
    ax[0][0].set_xticklabels("")
    # ax[0][0].set_ylabel(
    #     r"${\rm Re\,}\tilde{\varepsilon}_{xx}/{\rm Re\,}{\varepsilon^f_{xx}(0)}$"
    # )
    # ax[1][0].set_ylabel(r"$\tan \tilde{\delta}_{xx}/ \tan {\delta^f_{xx}} (0)$")
    #
    ax[0][0].set_ylabel("normalized permittivity")
    ax[1][0].set_ylabel("normalized loss tangent")

    colors = sns.color_palette()
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

    custom_lines0 = [
        Line2D([0], [0], color=colors[0]),
        Line2D([0], [0], color=colors[1]),
        Line2D([0], [0], color=colors[2]),
        Line2D([0], [0], color=colors[3]),
        Line2D([0], [0], color=colors[4]),
    ]
    custom_lines1 = [
        Line2D([0], [0], color="k", ls="--"),
        Line2D([0], [0], color="k", ls="-"),
    ]

    ax[0][0].legend(custom_lines1, [r"uncoupled", r"coupled"], loc=1, ncol=1)
    leg = [r"$f={:.1f}$".format(f) for f in Fincl[iplot]]
    ax[1][1].legend(custom_lines0, leg, loc=3)

    ax[0][1].plot(Ebias, n_norm_lin, ls="--")
    ax[0][1].set_prop_cycle(None)
    ax[0][1].plot(Ebias, n_norm_nl, ls="-")
    # ax[0][1].set_xlabel("electric field (MV/m)")
    ax[0][1].set_xticklabels("")
    # ax[0][1].set_ylabel(r"$\tilde{n}/n^{\rm f}$")
    ax[0][1].set_ylabel("normalized tunability")
    ax[1][1].plot(Ebias, Aniso_factor_lin, ls="--")
    ax[1][1].set_prop_cycle(None)
    ax[1][1].plot(Ebias, Aniso_factor_nl, ls="-")
    ax[1][1].set_xlabel("bias electric field (MV/m)")

    # ax[1][1].set_ylabel(
    #     r"$\tilde{\rho} = \tilde{\varepsilon}_{xx}/\tilde{\varepsilon}_{yy}$"
    # )

    ax[1][1].set_ylabel("anisotropy factor")

    ax[0][0].set_ylim((0.1, 0.9))
    subplot_id(ax=ax[0][0], id="a")
    ax[1][0].set_ylim((0.75, 1.01))
    subplot_id(ax=ax[1][0], id="b")
    ax[0][1].set_ylim((0.9, 1.41))
    subplot_id(ax=ax[0][1], id="c")
    ax[1][1].set_ylim((0.37, 1.03))
    subplot_id(ax=ax[1][1], id="d")
    # for a1 in [0,1]:
    #     for a2 in [0,1]:
    #         ax[a1][a2].set_xlim((0, 2))


def plot_eff_cqf(fig, ax):
    ax.plot(Ebias, Knorm_lin, ls="--")
    ax.plot(Ebias, K_bulk, "k-")
    ax.set_prop_cycle(None)
    ax.plot(Ebias, Knorm_nl, ls="-")
    ax.set_yscale("log")
    # ax.plot(Ebias, K_nl, ls="-")
    ax.set_xlabel("electric field (MV/m)")
    # plt.plot(Ebias, K_bulk, "--k")


if __name__ == "__main__":

    norm_eps, norm_eps_nl, norm_loss, norm_loss_nl, n_norm_lin, n_norm_nl, Aniso_factor_lin, Aniso_factor_nl = (
        load_results()
    )
    plt.close("all")
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(6.2, 4.2))
    plot_eff_par(fig, ax)
    plt.tight_layout()
    fig.savefig(os.path.join(rootdir, "data", "figures", "effective_params_per.eps"))
    #
    # fig, ax = plt.subplots()
    # plot_eff_cqf(fig, ax)
