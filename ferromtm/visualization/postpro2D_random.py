from ferromtm.visualization.plots import *
from ferromtm.models.coupled2D import *

n_bulk = tunability(epsilonr_ferroelectric(Ebias), epsilonr_ferroelectric(0))
K_bulk = cqf(epsilonr_ferroelectric(Ebias), epsilonr_ferroelectric(0))
eps_00 = epsilonr_ferroelectric(0).real
tandelta_00 = -epsilonr_ferroelectric(0).imag / epsilonr_ferroelectric(0).real


def load_arch(E_bias, f, isample, coupling=True):
    name = "rand_f_{:.2f}_E_{:.2f}_sample_{}".format(f, E_bias, isample)
    if not coupling:
        name += "_uncoupled"
    saveddir = os.path.join(data_folder, "rand_circ_rods")
    filename = os.path.join(saveddir, name + ".npz")
    arch = np.load(filename)
    eps_hom = arch["eps_hom"]
    epsi = arch["epsi"]
    E = arch["E"]
    return epsi, E, eps_hom


iplot = [0, 1, 2, 3, 4]
Nbsamples = 21


def load_results_():
    Eps_eff_lin = np.zeros((Nbsamples, nE, nF), dtype=complex)
    Eps_eff_nl = np.zeros((Nbsamples, nE, nF), dtype=complex)
    Aniso_factor_lin = np.zeros((Nbsamples, nE, nF), dtype=float)
    Aniso_factor_nl = np.zeros((Nbsamples, nE, nF), dtype=float)
    #
    sample_list = range(Nbsamples)

    for iF, f in enumerate(Fincl):
        # print(iF)
        for iE, E_bias in enumerate(Ebias):
            # print(iE)
            for isample in sample_list:
                # print(isample)
                print("E={:.2f}MV/m, f={:.2f}, sample {}".format(E_bias, f, isample))
                _, _, eps_hom = load_arch(E_bias, f, isample, coupling=True)
                # eigval, eigvec = np.linalg.eig(eps_hom)
                # eps_hom = np.diag(eigval)
                eps_eff_nl = eps_hom[0, 0]
                Eps_eff_nl[isample, iE, iF] = eps_eff_nl
                Aniso_factor_nl[isample, iE, iF] = aniso_factor(eps_hom)

                _, _, eps_hom_lin = load_arch(E_bias, f, isample, coupling=False)
                # eigval, eigvec = np.linalg.eig(eps_hom_lin)
                # eps_hom_lin = np.diag(eigval)
                eps_eff_lin = eps_hom_lin[0, 0]
                Eps_eff_lin[isample, iE, iF] = eps_eff_lin
                Aniso_factor_lin[isample, iE, iF] = aniso_factor(eps_hom_lin)

    Eps_eff_lin = Eps_eff_lin[:, :, iplot]
    Eps_eff_nl = Eps_eff_nl[:, :, iplot]
    Aniso_factor_lin = Aniso_factor_lin[:, :, iplot]
    Aniso_factor_nl = Aniso_factor_nl[:, :, iplot]

    norm_eps = Eps_eff_lin.real / eps_00
    norm_loss = -Eps_eff_lin.imag / Eps_eff_lin.real / tandelta_00
    norm_eps_nl = Eps_eff_nl.real / eps_00
    norm_loss_nl = -Eps_eff_nl.imag / Eps_eff_nl.real / tandelta_00

    n_lin = np.array([tunability(e, e[0, :]) for e in Eps_eff_lin])
    n_nl = np.array([tunability(e, e[0, :]) for e in Eps_eff_nl])
    K_lin = np.array([cqf(e, e[0, :]) for e in Eps_eff_lin])
    K_nl = np.array([cqf(e, e[0, :]) for e in Eps_eff_nl])

    n_norm_lin = np.zeros((Nbsamples, nE, nF))
    n_norm_nl = np.zeros((Nbsamples, nE, nF))
    # Knorm_lin = np.zeros((Nbsamples, nE, nF))
    # Knorm_nl = np.zeros((Nbsamples, nE, nF))
    for iF, f in enumerate(Fincl):
        for isample in sample_list:
            n_norm_lin[isample, :, iF] = n_lin[isample, :, iF] / n_bulk
            n_norm_nl[isample, :, iF] = n_nl[isample, :, iF] / n_bulk
    # Knorm_lin[isample, :, iF] = K_lin[isample, :, iF] / K_bulk
    # Knorm_nl[isample, :, iF] = K_nl[isample, :, iF] / K_bulk
    name = "random_case"
    saveddir = os.path.join(data_folder, "rand_circ_rods")
    filename = os.path.join(saveddir, name + ".npz")
    np.savez(
        filename,
        norm_eps=norm_eps,
        norm_eps_nl=norm_eps_nl,
        norm_loss=norm_loss,
        norm_loss_nl=norm_loss_nl,
        n_norm_lin=n_norm_lin,
        n_norm_nl=n_norm_nl,
        Aniso_factor_lin=Aniso_factor_lin,
        Aniso_factor_nl=Aniso_factor_nl,
    )

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


def load_results(retrieve=False):
    print("Loading results...")
    if retrieve:
        return load_results_()
    else:
        name = "random_case"
        saveddir = os.path.join(data_folder, "rand_circ_rods")
        filename = os.path.join(saveddir, name + ".npz")
        arch = np.load(filename)
        norm_eps = arch["norm_eps"]
        norm_eps_nl = arch["norm_eps_nl"]
        norm_loss = arch["norm_loss"]
        norm_loss_nl = arch["norm_loss_nl"]
        n_norm_lin = arch["n_norm_lin"]
        n_norm_nl = arch["n_norm_nl"]
        Aniso_factor_lin = arch["Aniso_factor_lin"]
        Aniso_factor_nl = arch["Aniso_factor_nl"]

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


def plot_eff_par(fig, ax, lin=False):
    colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    Ebias_plot = np.tile(Ebias, Nbsamples)

    if lin:
        data_1 = norm_eps
        data_2 = norm_loss
        data_3 = n_norm_lin
        data_4 = Aniso_factor_lin
    else:
        data_1 = norm_eps_nl
        data_2 = norm_loss_nl
        data_3 = n_norm_nl
        data_4 = Aniso_factor_nl
    for iF, f in enumerate(Fincl):
        clr = colors[iF]

        data = data_1[:, :, iF].ravel()
        sns.lineplot(x=Ebias_plot, y=data, ax=ax[0][0], ci="sd", color=clr)

        data = data_2[:, :, iF].ravel()
        sns.lineplot(x=Ebias_plot, y=data, ax=ax[1][0], ci="sd", color=clr)

        data = data_3[:, :, iF].ravel()
        sns.lineplot(x=Ebias_plot, y=data, ax=ax[0][1], ci="sd", color=clr)

        data = data_4[:, :, iF].ravel()
        sns.lineplot(x=Ebias_plot, y=data, ax=ax[1][1], ci="sd", color=clr)

    ax[1][0].set_xlabel("bias electric field (MV/m)")
    ax[0][0].set_xticklabels("")
    # ax[0][0].set_ylabel(
    #     r"${\rm Re\,}\tilde{\varepsilon}_{xx}/{\rm Re\,}{\varepsilon^f_{xx}(0)}$"
    # )
    # ax[1][0].set_ylabel(r"$\tan \tilde{\delta}_{xx}/ \tan {\delta^f_{xx}} (0)$")

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

    # ax[0][0].legend(custom_lines1, [r"uncoupled", r"coupled"], loc=1, ncol=1)
    leg = [r"$f={:.1f}$".format(f) for f in Fincl[iplot]]
    ax[0][1].legend(custom_lines0, leg, loc=3)
    ax[0][1].set_xticklabels("")
    # ax[0][1].set_ylabel(r"$\tilde{n}/n^{\rm f}$")

    ax[0][1].set_ylabel("normalized tunability")
    ax[1][1].set_xlabel("bias electric field (MV/m)")

    # ax[1][1].set_ylabel(
    #     r"$\tilde{\rho} = \tilde{\varepsilon}_{xx}/\tilde{\varepsilon}_{yy}$"
    # )

    ax[1][1].set_ylabel("anisotropy factor")

    ax[0][0].set_ylim((0.1, 0.9))
    subplot_id(ax=ax[0][0], id="a")
    ax[1][0].set_ylim((0.55, 1.01))
    subplot_id(ax=ax[1][0], id="b")
    ax[0][1].set_ylim((0.65, 1.1))
    subplot_id(ax=ax[0][1], id="c")
    ax[1][1].set_ylim((0.5, 1.14))
    subplot_id(ax=ax[1][1], id="d")
    #
    # for ax_ in ax.ravel():
    #     ax_.set_rasterized(True)


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

    norm_eps, norm_eps_nl, norm_loss, norm_loss_nl, n_norm_lin, n_norm_nl, Aniso_factor_lin, Aniso_factor_nl = load_results(
        retrieve=False
    )

    plt.close("all")
    fig, ax = plt.subplots(ncols=2, nrows=2)
    plot_eff_par(fig, ax, lin=True)

    fig, ax = plt.subplots(ncols=2, nrows=2)
    plot_eff_par(fig, ax, lin=False)

# data = norm_eps[:, :, iF].ravel()
# sns.lineplot(x=Ebias_plot, y=data, ax=ax[0][0], ci="sd", color=clr)
# ax[0, 0].lines[2 * iF].set_linestyle("--")
#
# data = norm_loss[:, :, iF].ravel()
# sns.lineplot(x=Ebias_plot, y=data, ax=ax[1][0], ci="sd", color=clr)
# ax[1, 0].lines[2 * iF].set_linestyle("--")
#
# data = n_norm_lin[:, :, iF].ravel()
# sns.lineplot(x=Ebias_plot, y=data, ax=ax[0][1], ci="sd", color=clr)
# ax[0, 1].lines[2 * iF].set_linestyle("--")
#
# data = Aniso_factor_lin[:, :, iF].ravel()
# sns.lineplot(x=Ebias_plot, y=data, ax=ax[1][1], ci="sd", color=clr)
#
# ax[1, 1].lines[2 * iF].set_linestyle("--")


# plot_eff_par(fig, ax)
#
# fig, ax = plt.subplots()
# plot_eff_cqf(fig, ax)
