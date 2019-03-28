from ferromtm.models.bst import *
from ferromtm.visualization.plots import *


if __name__ == "__main__":
    fit_params = np.load(os.path.join(rootdir, "data", "fit_params.npz"))["fit_params"]
    meas = np.load(os.path.join(rootdir, "data", "measurements.npz"))
    E_exp = meas["E_exp"]
    eps_norm_exp = meas["eps_norm_exp"]
    eps_exp_0 = meas["eps_exp_0"]

    E_fit = np.linspace(-6, 6, 101)
    eps_norm_fit = []
    for i in range(2):
        alpha_opt = fit_params[i]
        epsf0 = eps_exp_0[i]
        eps_fit = eps_norma_model(E_fit, *alpha_opt)
        eps_norm_fit.append(eps_fit.real)

    plt.close("all")
    fig = plt.figure(figsize=(5, 3))
    plt.plot(E_exp, eps_norm_exp[0], "o", label="measured, static", color="#3e8b58")
    plt.plot(E_fit, eps_norm_fit[0], "--", color="#3e8b58", label="fit, static")
    plt.plot(
        E_exp, eps_norm_exp[1], "o", label="measured, $f=3.8$ GHz", color="#eb931e"
    )
    plt.plot(E_fit, eps_norm_fit[1], "--", color="#eb931e", label="fit, $f=3.8$ GHz")
    plt.ylabel(r"normalized permittivity")
    # $\varepsilon^{\rm f}/\varepsilon^{\rm f}(E=0)$
    plt.xlabel("electric field (MV/m)")
    plt.legend()

    plt.tight_layout()
    plt.ylim((0.1, 1.04))
    fig.savefig(os.path.join(rootdir, "data", "figures", "epsilon_fit.eps"))
