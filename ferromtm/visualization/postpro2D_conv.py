from ferromtm.models.coupled2D import *
from ferromtm.visualization.plots import *
from matplotlib.ticker import MaxNLocator


n_bulk = tunability(epsilonr_ferroelectric(Ebias), epsilonr_ferroelectric(0))
K_bulk = cqf(epsilonr_ferroelectric(Ebias), epsilonr_ferroelectric(0))
eps_00 = epsilonr_ferroelectric(0).real
tandelta_00 = -epsilonr_ferroelectric(0).imag / epsilonr_ferroelectric(0).real

cv_dir_ = "circ_rods"
cv_dir_ = "rand_circ_rods"
cv_dir = os.path.join(data_folder, cv_dir_, "convergence")


def load_arch(iter):
    fname = "cv_iter_{}.npz".format(iter)
    filename = os.path.join(cv_dir, fname)
    arch = np.load(filename)
    eps_hom = arch["eps_hom"]
    epsi = arch["epsi_map"]
    E = arch["E_map"]
    return epsi, E, eps_hom


ncv = 13  # periodic
ncv = 15  # random


def load_results():
    Eps_eff_xx = np.zeros(ncv, dtype=complex)
    Eps_eff_yy = np.zeros(ncv, dtype=complex)
    Aniso_factor = np.zeros(ncv, dtype=float)
    E_cv = []
    epsi_cv = []
    for icv in range(ncv):
        epsi, E, eps_hom = load_arch(icv)
        Eps_eff_xx[icv] = eps_hom[0, 0]
        Eps_eff_yy[icv] = eps_hom[1, 1]
        Aniso_factor[icv] = aniso_factor(eps_hom)
        E_cv.append(E)
        epsi_cv.append(epsi)

    norm_eps_xx = Eps_eff_xx.real / eps_00
    norm_loss_xx = -Eps_eff_xx.imag / Eps_eff_xx.real / tandelta_00
    norm_eps_yy = Eps_eff_yy.real / eps_00
    norm_loss_yy = -Eps_eff_yy.imag / Eps_eff_yy.real / tandelta_00
    return (norm_eps_xx, norm_loss_xx, norm_eps_yy, norm_loss_yy, E_cv, epsi_cv)


E0x, E0y = 2, 0
normE0 = 2


def plotEarrows(ax, Enewx, Enewy):
    Enewx += E0x
    Enewy += E0y
    dspl = 2 ** 0
    nx, ny = Enewx.shape
    dd = np.linspace(0, nx - 1, nx)
    X0, Y0 = np.meshgrid(dd, dd)
    X, Y = X0[::dspl, ::dspl], Y0[::dspl, ::dspl]
    X, Y = X0, Y0

    U0 = Enewx[::dspl, ::dspl].T
    V0 = Enewy[::dspl, ::dspl].T
    normarrow = np.sqrt(U0 ** 2 + V0 ** 2)
    U = U0
    V = V0
    ax.streamplot(X, Y, U, V, density=1, linewidth=1, color="#ababab")


def plot_E_map_conv(axtmp, i, label=None):
    epsp = epsi_cv[i][0].real
    epsp = epsp.T
    E = E_cv[i]
    Enewx, Enewy = E
    Enewx += E0x
    Enewy += E0y
    Enp = (Enewx ** 2 + Enewy ** 2) ** 0.5
    Enp = Enp.T / normE0
    Enp[epsp <= 3] = np.NaN
    plotEarrows(axtmp, Enewx, Enewy)
    Eax = axtmp.imshow(np.log10(Enp), cmap="YlOrRd", vmin=-1, vmax=1)
    axtmp.set_aspect("equal")
    axtmp.set_axis_off()
    cb = plt.colorbar(Eax, ax=axtmp)
    axtmp.set_title(r"electric field $E/E_0$, $i={0}$".format(i))

    return axtmp, cb


def plot_eps_map_conv(ax, i, label=None):

    epsp_ = epsi_cv[i][0].real

    epsp_ = epsp_.T
    epsp_[epsp_ <= 3] = np.NaN
    epsax = ax.imshow(epsp_, cmap="summer", vmax=eps_00.real)
    ax.set_aspect("equal")
    ax.set_axis_off()
    cb = plt.colorbar(epsax, ax=ax)
    ax.set_title(r"permittivity $\varepsilon_{{xx}}$, $i={0}$".format(i))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    return ax, cb


def plot_eps_conv(ax):

    ax.plot(
        norm_eps_xx,
        "-",
        color="#e79037",
        label=r"$xx$",
        # label=r"$\tilde{\varepsilon}_{xx}/\varepsilon^{\rm f}_{xx}(0)$",
    )
    ax.plot(
        norm_eps_yy,
        "--",
        color="#67a671",
        label=r"$yy$",
        # label=r"$\tilde{\varepsilon}_{yy}/\varepsilon^{\rm f}_{yy}(0)$",
    )
    ax.legend()
    ax.set_xticklabels("")
    ax.set_title("normalized permittivity")
    ax.set_ylim((0.14, 0.22))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    return ax


def plot_tand_conv(axtmp):
    axtmp.plot(norm_loss_xx, "-", color="#e79037", label=r"$xx$")
    # label=r"$\tan\tilde{\delta}_{xx}/\tan\delta^{\rm f}_{xx}(0)$")
    axtmp.plot(norm_loss_yy, "--", color="#67a671", label=r"$yy$")
    # label=r"$\tan\tilde{\delta}_{yy}/\tan\delta^{\rm f}_{yy}(0)$")
    axtmp.legend()
    axtmp.set_xlabel("iteration $i$")
    axtmp.set_title("normalized loss tangent")
    # axtmp.set_ylim((0.8, 0.97))
    axtmp.set_ylim((0.65, 0.85))
    return axtmp


if __name__ == "__main__":

    norm_eps_xx, norm_loss_xx, norm_eps_yy, norm_loss_yy, E_cv, epsi_cv = load_results()
    plt.close("all")
    fig, ax = plt.subplots(2, 3)
    plot_eps_conv(ax[0, 0])
    subplot_id(ax=ax[0, 0], id="a")
    plot_tand_conv(ax[1, 0])
    subplot_id(ax=ax[1, 0], id="b")
    i = 1
    plot_E_map_conv(ax[0, 1], i)
    subplot_id(ax=ax[0, 1], id="c")
    plot_eps_map_conv(ax[0, 2], i)
    subplot_id(ax=ax[0, 2], id="e")

    i = ncv - 1
    # i = 1 # random
    plot_E_map_conv(ax[1, 1], i)
    subplot_id(ax=ax[1, 1], id="d")
    plot_eps_map_conv(ax[1, 2], i)
    subplot_id(ax=ax[1, 2], id="f")
