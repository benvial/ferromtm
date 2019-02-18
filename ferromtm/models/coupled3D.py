#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


# import matplotlib.pyplot as plt
# plt.ion()
import numpy as np
from aotomat.tools.plottools import *
from pytheas.homogenization import *
from aotomat.electrostatics.per3D import femmodel as model_es
from pytheas.homogenization.twoscale3D import femmodel as model_hom
import importlib

importlib.reload(model_es)
importlib.reload(model_hom)


from ferromtm.tools.utils import *

pi = np.pi

eps_incl = 3


def compute_hom_pb(fem_hom, epsi):
    """Computes the homogenization problem.

    Parameters
    ----------
    fem_hom : instance of pytheas.homogenization.twoscale3D.femmodel.TwoScale3D class
        The two-scale FEM problem
    epsi : array-like
        The permittivity tensor diagonal elements, (xx, yy and zz components).
        Size is (3*nvar) where nvar is the number of points where the permittivity
        is defined in the FEM model (ie len(fem_hom.des[0])).

    Returns
    -------
    tuple (eps_hom, fem_hom)
        eps_hom: the permittivity tensor (3*3 numpy array)
        fem_hom: instance of pytheas.homogenization.twoscale3D.femmodel.TwoScale3D class

    """

    make_pos_tensor_eps(fem_hom, epsi)
    fem_hom.compute_solution()
    eps_hom = fem_hom.compute_epsilon_eff()
    return eps_hom, fem_hom


def compute_elstat_pb(fem_es, epsi):
    make_pos_tensor_eps(fem_es, epsi)
    fem_es.compute_solution()
    E = fem_es.postpro_electrostatic_field()
    return E, fem_es


def coupling_loop(fem_es, epsi, tol=1e-2, verbose=True):
    conv = False
    Eold = fem_es.E_static
    while not conv:
        E, fem_es = compute_elstat_pb(fem_es, epsi)
        epsi = epsilonr_ferroelectric(E.real)
        cv = np.mean(np.abs(E - Eold) ** 2) ** 0.5
        conv = cv < tol
        Eold = E
        if verbose:
            print("error: ", cv)
    return E, epsi, fem_es


def init_es(f, E_bias):
    r = (3 * f / (4 * pi)) ** (1 / 3)
    #####################################
    # # electrostatics
    #####################################
    fem_es = model_es.FemModel()
    fem_es.gmsh_verbose = 0  #: str: Gmsh verbose (int between 0 and 4)
    fem_es.getdp_verbose = 0  #: str: GetDP verbose (int between 0 and 4)
    fem_es.python_verbose = 0  #: str: GetDP verbose (int between 0 and 1)
    fem_es.parmesh = 4
    fem_es.parmesh_incl = 4
    fem_es.type_des = "elements"
    fem_es.inclusion_flag = True
    fem_es.dx = 1
    fem_es.dy = 1
    fem_es.dz = 1
    ## principal axes of the ellipsoidal inclusion
    fem_es.ax = r
    fem_es.ay = r
    fem_es.az = r
    fem_es.eps_incl = eps_incl
    fem_es.eps_host = eps_incl
    fem_es.E_static = E_bias
    fem_es.coupling_flag = True
    fem_es.tmp_dir = "./estat"
    fem_es.rm_tmp_dir()
    fem_es.initialize()
    fem_es.make_mesh()
    return fem_es


def init_hom(fem_es):
    #####################################
    # # homogenization
    #####################################
    fem_hom = model_hom.TwoScale3D()
    fem_hom.gmsh_verbose = fem_es.gmsh_verbose
    fem_hom.getdp_verbose = fem_es.getdp_verbose
    fem_hom.python_verbose = fem_es.python_verbose
    fem_hom.parmesh = fem_es.parmesh
    fem_hom.parmesh_incl = fem_es.parmesh_incl
    fem_hom.type_des = fem_es.type_des
    fem_hom.inclusion_flag = fem_es.inclusion_flag
    fem_hom.dx = fem_es.dx
    fem_hom.dy = fem_es.dy
    fem_hom.dz = fem_es.dz
    ## principal axes of the ellipsoidal inclusion
    fem_hom.ax = fem_es.ax
    fem_hom.ay = fem_es.ay
    fem_hom.az = fem_es.az
    fem_hom.eps_incl = fem_es.eps_incl
    fem_hom.eps_host = fem_es.eps_host
    fem_hom.tmp_dir = "./hom"
    fem_hom.coupling_flag = fem_es.coupling_flag
    fem_hom.rm_tmp_dir()

    fem_hom.initialize()
    fem_hom.make_mesh()
    return fem_hom


def main(f, E_bias, coupling=True):
    fem_es = init_es(f, E_bias)
    nvar = len(fem_es.des[0])
    eps_f = epsilonr_ferroelectric(E_bias)
    eps_f0 = epsilonr_ferroelectric(0)
    id = np.ones(nvar)
    epsi = id * eps_f, id * eps_f0, id * eps_f0
    E = E_bias * id
    if E_bias != 0:
        if coupling:
            E, epsi, fem_es = coupling_loop(fem_es, epsi)
    fem_hom = init_hom(fem_es)
    eps_hom, fem_hom = compute_hom_pb(fem_hom, epsi)
    return eps_hom, epsi, E, fem_hom, fem_es


def efield_mean_mg(E0, f, epsi, epsh):
    """Mean E field in the background for an isolated sphere"""
    return E0 / (1 - f + 3 * f * epsh / (epsi + 2 * epsh))


def coupling_loop_analytical(E, f, eps_incl, eps_f, verbose=False):
    """
    mean_field(E, epsilon)  <---> epsilon(E)
    """
    conv = False
    E0 = E
    Eold = E
    i = 0
    while not conv:
        E = efield_mean_mg(E0, f, eps_incl, eps_f)
        eps_f = epsilonr_ferroelectric(E.real)
        cv = np.mean(np.abs(E - Eold) ** 2) ** 0.5
        if verbose:
            print("error: ", cv)
        conv = cv < 1e-2
        Eold = E
        i += 1
        if i > 100:
            break
    return E, eps_f


if __name__ == "__main__":
    fmax = 4 / 3 * pi * 0.5 ** 3
    f = 0.4
    nbias = 10
    bias = np.linspace(0, 2, nbias)

    eps_f = epsilonr_ferroelectric(bias)
    eps_f0 = epsilonr_ferroelectric(0)
    eps_harm = harmonic_mean(f, eps_incl, eps_f)
    eps_ar = arithmetic_mean(f, eps_incl, eps_f)

    eps_c, eps_s = [], []
    epsilon_c, epsilon_s = [], []
    E_c, E_s = [], []
    for E_bias in bias:
        print("")
        print("*" * 33)
        print("Bias = {} MV/m".format(E_bias))
        print("*" * 33)
        print("")
        print("Coupled --------------")
        eps_hom, epsi, E, fem_hom, fem_es = main(f, E_bias, coupling=True)
        print("eps_hom")
        matprint(eps_hom)
        print("Uncoupled ------------")
        eps_hom_uc, epsi_uc, E_uc, fem_hom, fem_es = main(f, E_bias, coupling=False)
        print("eps_hom_uc")
        matprint(eps_hom_uc)
        eps_c.append(eps_hom)
        eps_s.append(eps_hom_uc)
        epsilon_c.append(epsi)
        epsilon_s.append(epsi_uc)
        E_c.append(E)
        E_s.append(E_uc)

    eps_bias_c, eps_bias_s = (
        np.zeros((3, nbias), dtype=complex),
        np.zeros((3, nbias), dtype=complex),
    )
    eps_bias_c_mean, E_bias_c_mean = (
        np.zeros((3, nbias), dtype=complex),
        np.zeros((3, nbias), dtype=complex),
    )
    for i in range(3):
        eps_bias_c[i, :] = np.array([_[i, i] for _ in eps_c])
        eps_bias_s[i, :] = np.array([_[i, i] for _ in eps_s])
        eps_bias_c_mean[i, :] = np.array([np.mean(_[i]) for _ in epsilon_c])
        # E_bias_c_mean[i,:] = np.array([np.mean(_[i]) for _ in E_c])

    eps_mg = [maxwell_garnett(f, eps_incl, eps_f)]
    eps_mg.append(np.ones(nbias) * eps_mg[0][0])
    eps_mg_mean = maxwell_garnett(f, eps_incl, eps_bias_c_mean)

    eps_ana_c, E_ana_c = [], []
    for E in bias:
        eps_ = epsilonr_ferroelectric(E)
        E, eps_ = coupling_loop_analytical(E, f, eps_incl, eps_)
        eps_ana_c.append(eps_)
        E_ana_c.append(E)
    eps_ana_c = np.array(eps_ana_c)
    eps_mg_ana_c = [maxwell_garnett(f, eps_incl, eps_ana_c)]
    eps_mg_ana_c.append(np.ones(nbias) * eps_mg_ana_c[0][0])

    plt.close("all")
    fig, ax = plt.subplots(2, figsize=(4, 8))
    ax[0].plot(bias, eps_bias_c[0], "-r", label="coupled $\parallel$")
    ax[0].plot(bias, eps_bias_s[0], "-b", label="uncoupled $\parallel$")
    # ax[0].plot(bias, eps_mg[0], "-g", label="Maxwell-Garnett $\parallel$")
    # ax[0].plot(bias, eps_mg_mean[0], "-m", label="numerical mean epsi MG $\parallel$")
    # ax[0].plot(bias, eps_mg_ana_c[0], "-c", label="analytical mean epsi MG $\parallel$")

    ax[0].legend()
    ax[1].plot(bias, eps_bias_c[1], "-r", label="coupled $\perp$")
    ax[1].plot(bias, eps_bias_s[1], "-b", label="uncoupled $\perp$")
    # ax[1].plot(bias, eps_mg[1], "-g", label="Maxwell-Garnett $\perp$")
    # ax[1].plot(bias, eps_mg_mean[1], "-m", label="numerical mean epsi MG $\perp$")
    # ax[1].plot(bias, eps_mg_ana_c[1], "-c", label="analytical mean epsi MG $\perp$")

    ax[1].legend()
    ax[1].set_xlabel("Bias field $E$ (MV/m)")
    ax[0].set_title("effective permittivity, $f = {}$".format(f))

    ## loss tangent
    tand_bias_c = losstan(eps_bias_c)
    tand_bias_s = losstan(eps_bias_s)
    tand_mg = losstan(eps_mg)
    tand_mg_mean = losstan(eps_mg_mean)
    tand_mg_ana_c = losstan(eps_mg_ana_c)

    fig, ax = plt.subplots(2, figsize=(4, 8))
    ax[0].plot(bias, tand_bias_c[0], "-r", label="coupled $\parallel$")
    ax[0].plot(bias, tand_bias_s[0], "-b", label="uncoupled $\parallel$")
    # ax[0].plot(bias, tand_mg[0], "-g", label="Maxwell-Garnett $\parallel$")
    # ax[0].plot(bias, tand_mg_mean[0], "-m", label="numerical mean epsi MG $\parallel$")
    # ax[0].plot(bias, tand_mg_ana_c[0], "-c", label="analytical mean epsi MG $\perp$")

    ax[0].legend()
    ax[1].plot(bias, tand_bias_c[1], "-r", label="coupled $\perp$")
    ax[1].plot(bias, tand_bias_s[1], "-b", label="uncoupled $\perp$")
    # ax[1].plot(bias, tand_mg[1], "-g", label="Maxwell-Garnett $\perp$")
    # ax[1].plot(bias, tand_mg_mean[1], "-m", label="numerical mean epsi MG $\perp$")
    # ax[1].plot(bias, tand_mg_ana_c[1], "-c", label="analytical mean epsi MG $\perp$")

    ax[1].legend()
    ax[1].set_xlabel("Bias field $E$ (MV/m)")
    ax[0].set_title("effective loss tangent, $f = {}$".format(f))

    ## tunabilities

    tunability_c = []
    for e in eps_bias_c:
        tunability_c.append(tunability(e, e[0]))
    tunability_s = []
    for e in eps_bias_s:
        tunability_s.append(tunability(e, e[0]))
    tunability_mg = []
    for e in eps_mg:
        tunability_mg.append(tunability(e, e[0]))
    tunability_mg_mean = []
    for e in eps_mg_mean:
        tunability_mg_mean.append(tunability(e, e[0]))
    tunability_mg_ana_c = []
    for e in eps_mg_ana_c:
        tunability_mg_ana_c.append(tunability(e, e[0]))

    tunability_bulk = [tunability(eps_f, eps_f0)]
    tbx = [1 for _ in range(nbias)]
    tunability_bulk.append(tbx)
    tunability_bulk.append(tbx)

    fig, ax = plt.subplots(2, figsize=(4, 8))
    ax[0].plot(bias, tunability_c[0], "-r", label="coupled $\parallel$")
    ax[0].plot(bias, tunability_s[0], "-b", label="uncoupled $\parallel$")
    # ax[0].plot(bias, tunability_mg[0], "-g", label="Maxwell-Garnett $\parallel$")
    # ax[0].plot(
    #     bias, tunability_mg_mean[0], "-m", label="numerical mean epsi MG $\parallel$"
    # )
    # ax[0].plot(
    #     bias, tunability_mg_ana_c[0], "-c", label="analytical mean epsi MG $\parallel$"
    # )

    ax[0].plot(bias, tunability_bulk[0], "-k", label="bulk $\parallel$")
    ax[0].legend()
    ax[1].plot(bias, tunability_c[1], "-r", label="coupled $\perp$")
    ax[1].plot(bias, tunability_s[1], "-b", label="uncoupled $\perp$")
    # ax[1].plot(bias, tunability_mg[1], "-g", label="Maxwell-Garnett $\perp$")
    # ax[1].plot(
    #     bias, tunability_mg_mean[1], "-m", label="numerical mean epsi MG $\perp$"
    # )
    # ax[1].plot(
    #     bias, tunability_mg_ana_c[1], "-c", label="analytical mean epsi MG $\parallel$"
    # )

    ax[1].plot(bias, tunability_bulk[1], "-k", label="bulk $\perp$")
    ax[1].legend()
    ax[1].set_xlabel("Bias field $E$ (MV/m)")
    ax[0].set_title("tunability, $f = {}$".format(f))
    #

    # commutation quality factor

    K_c = []
    for e in eps_bias_c:
        K_c.append(cqf(e, e[0]))
    K_s = []
    for e in eps_bias_s:
        K_s.append(cqf(e, e[0]))

    K_mg = []
    for e in eps_mg:
        K_mg.append(cqf(e, e[0]))
    K_mg_mean = []
    for e in eps_mg_mean:
        K_mg_mean.append(cqf(e, e[0]))
    K_mg_ana_c = []
    for e in eps_mg_ana_c:
        K_mg_ana_c.append(cqf(e, e[0]))

    K_bulk = [cqf(eps_f, eps_f0)]
    Ktmp = [0 for _ in range(nbias)]
    K_bulk.append(Ktmp)
    K_bulk.append(Ktmp)

    fig, ax = plt.subplots(2, figsize=(4, 8))
    ax[0].plot(bias, K_c[0], "-r", label="coupled $\parallel$")
    ax[0].plot(bias, K_s[0], "-b", label="uncoupled $\parallel$")
    # ax[0].plot(bias, K_mg[0], "-g", label="Maxwell-Garnett $\parallel$")
    # ax[0].plot(
    #     bias, K_mg_mean[0], "-m", label="numerical mean epsi MG $\parallel$"
    # )
    # ax[0].plot(
    #     bias, K_mg_ana_c[0], "-c", label="analytical mean epsi MG $\parallel$"
    # )
    ax[0].plot(bias, K_bulk[0], "-k", label="bulk $\parallel$")
    ax[0].legend()
    ax[1].plot(bias, K_c[1], "-r", label="coupled $\perp$")
    ax[1].plot(bias, K_s[1], "-b", label="uncoupled $\perp$")
    # ax[1].plot(bias, K_mg[1], "-g", label="Maxwell-Garnett $\perp$")
    # ax[1].plot(
    #     bias, K_mg_mean[1], "-m", label="numerical mean epsi MG $\perp$"
    # )
    # ax[1].plot(
    #     bias, K_mg_ana_c[1], "-c", label="analytical mean epsi MG $\parallel$"
    # )

    ax[1].plot(bias, K_bulk[1], "-k", label="bulk $\perp$")
    ax[1].legend()
    ax[1].set_xlabel("Bias field $E$ (MV/m)")
    ax[0].set_title("commutation quality factor, $f = {}$".format(f))
