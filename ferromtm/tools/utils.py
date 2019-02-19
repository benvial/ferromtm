#!/usr/bin/env python

import numpy as np


def make_pos_tensor_eps(fem, epsi, interp=False):
    fem.path_pos = ""
    for epsi_comp, comp in zip(epsi, ["xx", "yy", "zz"]):
        if interp:
            feps = fem.make_fdens(epsi_comp)
            eps_des = feps(fem.des[1])
        else:
            eps_des = epsi_comp
        fem.path_pos += " " + fem.make_eps_pos(
            fem.des[0], eps_des, posname="eps_des_" + comp
        )


def tunability(eps_f, eps_f0):
    return eps_f0.real / eps_f.real


def losstan(epsi):
    return np.abs(np.imag(epsi) / np.real(epsi))


def cqf(eps, eps0):
    n = tunability(eps, eps0)
    td = eps.imag / eps.real
    td0 = eps0.imag / eps0.real
    return (n - 1) ** 2 / (n * td0 * td)


def aniso_factor(epsilon):
    return epsilon[0, 0].real / epsilon[1, 1].real


def matprint(mat):
    print("")
    for x in mat:
        print("{:.3f}    {:.3f}      {:.3f}".format(x[0], x[1], x[2]))
    print("")


def subplot_id(id="a", ax=None):
    ax.text(
        -0.25, 1.05, id, transform=ax.transAxes, fontsize="large", fontweight="bold"
    )
