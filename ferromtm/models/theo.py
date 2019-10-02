#!/usr/bin/env python

import numpy as np

from ferromtm import rootdir
from ferromtm.models.bst import *

epsilon0 = 8.8541878128e-12

# fname = "drilled vs undrilled"
from ferromtm import rootdir

fname = os.path.join(
    rootdir, "data/theo/drilled bst/drilled vs undrilled simalr feild pe"
)
f_ = np.loadtxt(fname)

ind0 = 180
ind1 = 311
ind2 = 465
ind3 = 601
#
# ind0 = 0
# ind1 = 601

# f = np.vstack((f_[ind0:ind1, :], np.flipud(f_[ind2:ind3, :])))
# f = np.vstack((f_[ind0:ind1],(f_[ind2:ind3])))

# f = f_[ind0:ind1, :]
# f = np.flipud(f_[ind2:ind3, :])
f = f_

Eh = f[:, 0]
Ph = f[:, 1]
E = f[:, 2]
P = f[:, 3]

# ind0 = 0
# ind1 = 311
# #
# # ind0 = 0
# # ind1 = 601
#
# Eh = f[ind0:ind1,0]
# Ph = f[ind0:ind1,1]
# E = f[ind0:ind1,2]
# P = f[ind0:ind1,3]
#
i0 = 301
eps = np.gradient(P) / np.gradient(E) * 1e-6 / epsilon0
eps0 = eps[i0]
eps = np.hstack((eps[ind0:ind1], np.flipud(eps[ind2:ind3])))
E_ = np.hstack((E[ind0:ind1], np.flipud(E[ind2:ind3])))
# eps0 = np.max(eps)
# eps0 = eps[601-ind0]
eps /= eps0
epsh = np.gradient(Ph) / np.gradient(Eh) * 1e-6 / epsilon0
epsh0 = epsh[i0]
epsh = np.hstack((epsh[ind0:ind1], np.flipud(epsh[ind2:ind3])))
Eh_ = np.hstack((Eh[ind0:ind1], np.flipud(Eh[ind2:ind3])))
epsh /= epsh0


E_Theo = E_
E_Theo_h = Eh_
eps_Theo = eps
eps_Theo_h = epsh
eps_Theo_0 = eps0
eps_Theo_h_0 = epsh0


def fit(filename, fitname, par0=(1, 1)):
    meas = np.load(os.path.join(rootdir, "data", filename))
    E_exp = meas["E_exp"]
    eps_norm_exp = meas["eps_norm_exp"]
    eps_exp_0 = meas["eps_exp_0"]
    fit_params = []
    alpha_opt = retrieve_params(E_exp, eps_norm_exp, par0=par0)
    fit_params = list(alpha_opt)
    np.savez(os.path.join(rootdir, "data", fitname), fit_params=fit_params)
    fit_params.append(eps_exp_0)
    print(fit_params)
    return fit_params


measurements_name = "measurements_theo.npz"
np.savez(
    os.path.join(rootdir, "data", measurements_name),
    eps_norm_exp=eps,
    E_exp=E_,
    eps_exp_0=eps0,
)

fit_params = fit("measurements_theo.npz", "fit_params_theo.npz", par0=(0.1, 0))


def epsf_real(E_applied, fit_params):
    alpha, beta, eps00 = fit_params
    return eps_norma_model(E_applied, alpha, beta) * eps00


E_applied_i = np.linspace(-3, 3, 1001)
epsilonr_ferroelectric_i = epsf_real(E_applied_i, fit_params)


def epsilonr_ferroelectric(E_applied, tandelta=1e-2):
    epsi = epsilonr_ferroelectric_i
    return np.interp(E_applied, E_applied_i, epsi) * (1 - tandelta * 1j)


if __name__ == "__main__":

    from aotomat.tools.plottools import *

    col1 = "#6078cf"
    col2 = "#c13f3f"
    plt.close("all")

    plt.figure()
    plt.plot(E, P, "o", label="bulk", color=col2, alpha=0.3)
    plt.plot(Eh, Ph, "o", label="mtm", color=col1, alpha=0.3)
    plt.xlabel("$E$ (kV/mm)")
    plt.ylabel("$P$ (kV/mm)")
    plt.title("PE loop")
    plt.legend()
    # plt.figure()
    # plt.plot(E, P,"-", label="bulk")
    # plt.plot(Eh,Ph,"--", label="mtm")
    # plt.xlabel("$E$ (kV/mm)")
    # plt.ylabel("$P$")
    # plt.legend()
    plt.figure()
    plt.plot(E_, eps, "o", label="bulk", color=col2, alpha=0.3)
    plt.plot(Eh_, epsh, "o", label="mtm", color=col1, alpha=0.3)
    plt.xlabel("$E$ (kV/mm)")
    plt.ylabel("normalized permittivity")

    plt.plot(
        E_applied_i,
        epsilonr_ferroelectric_i / epsilonr_ferroelectric(0),
        "--",
        label="fit",
        color=col2,
    )
    plt.legend()
