import numpy as np
import os
import scipy as sc
from scipy.optimize import minimize
from ferromtm import rootdir


def eps_norma_model(E0, alpha, beta):
    E0_ = np.asarray(E0)
    E_ = E0_.ravel()
    out = []
    for E in E_:
        poly = [alpha / 3, 0, 1, -E]
        poly = [beta / 5, 0, alpha / 3, 0, 1, -E]
        x = np.roots(poly)
        P0 = x[np.isreal(x)][0].real
        D = 1 + alpha * P0 ** 2
        D = 1 + alpha * P0 ** 2 + beta * P0 ** 4
        out.append(1 / D * (1 - 0 * 1j))
    return np.asarray(out).reshape(E0_.shape)


# def fun_fit(alpha, E, eps_meas_n):
#     eps_model = eps_norma_model(E, alpha)
#     err = np.abs(eps_model - eps_meas_n) ** 2
#     return np.mean(err) ** 0.5


def fun_fit(par, E, eps_meas_n):
    alpha, beta = par
    eps_model = eps_norma_model(E, alpha, beta)
    err = np.abs(eps_model - eps_meas_n) ** 2
    return np.mean(err) ** 0.5


def retrieve_params(E, eps_meas_n, par0=(1, 1)):
    cons = {"type": "ineq", "fun": lambda x: np.array([x[0]])}
    opt = minimize(
        fun_fit, par0, args=(E, eps_meas_n), options={"disp": True}, constraints=cons
    )
    xopt = opt["x"]
    print("xopt = ", xopt)
    return xopt


def fit(filename, fitname):
    meas = np.load(os.path.join(rootdir, "data", filename))
    E_exp = meas["E_exp"]
    eps_norm_exp = meas["eps_norm_exp"]
    eps_exp_0 = meas["eps_exp_0"]
    fit_params = []
    for i in range(2):
        alpha_opt = retrieve_params(E_exp, eps_norm_exp[i])
        fit_params.append(alpha_opt)
    np.savez(os.path.join(rootdir, "data", fitname), fit_params=fit_params)
    return fit_params


def epsf_real(E_applied, dc=False):
    if dc:
        eps00 = 3050
        alpha, beta = 0.11957429, 0.02415848
    else:
        eps00 = 165
        alpha, beta = 0.2403613, 0.07910162
    return eps_norma_model(E_applied, alpha, beta) * eps00


E_applied_i = np.linspace(-10, 10, 1001)
epsilonr_ferroelectric_i = epsf_real(E_applied_i)
epsilonr_ferroelectric_i_dc = epsf_real(E_applied_i, dc=True)


def epsilonr_ferroelectric(E_applied, tandelta=1e-2, dc=False):
    if dc:
        epsi = epsilonr_ferroelectric_i_dc
    else:
        epsi = epsilonr_ferroelectric_i
    return np.interp(E_applied, E_applied_i, epsi) * (1 - tandelta * 1j)


if __name__ == "__main__":

    fit("measurements.npz", "fit_params.npz")
