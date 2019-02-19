#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT

from ferromtm import __file__

import numpy as np

# from aotomat.tools.plottools import *
from pytheas.homogenization import *
from aotomat.electrostatics.per2D import femmodel as model_es
from pytheas.homogenization.twoscale2D import femmodel as model_hom
import importlib
from ferromtm.tools.utils import *
from pytheas.material import genmat
import os
import tempfile
from ferromtm.models.bst import epsilonr_ferroelectric
from aotomat.tools.parallelize import *

rootdir = os.path.dirname(os.path.dirname(__file__))
#
#
#
# from aotomat.material.bst import epsi_bst_norma
#
#
# def epsilonr_ferroelectric_(E_applied, tandelta=1e-2, eps00=120, sample=4):
#     return epsi_bst_norma(E_applied, sample) * eps00 * (1 - tandelta * 1j)
#
#
# E_applied_i = np.linspace(-10, 10, 1001)
# epsilonr_ferroelectric_i = epsilonr_ferroelectric_(E_applied_i)
#
# def epsilonr_ferroelectric(E_applied, tandelta=1e-2, eps00=120, sample=4, dc=False):
#     return np.interp(E_applied, E_applied_i, epsilonr_ferroelectric_i)
#
#


# importlib.reload(model_es)
# importlib.reload(model_hom)

pi = np.pi
eps_incl = 3
data_folder = os.path.join(rootdir, "data", "results")
mat_folder = "../data/mat/random/"


def ellipse(Rinclx, Rincly, rot_incl, x0, y0):
    c, s = np.cos(rot_incl), np.sin(rot_incl)
    Rot = np.array([[c, -s], [s, c]])
    nt = 360
    theta = np.linspace(-pi, pi, nt)
    x = Rinclx * np.sin(theta)
    y = Rincly * np.cos(theta)
    x, y = np.linalg.linalg.dot(Rot, np.array([x, y]))
    points = x + x0, y + y0
    return points


def circle(R, x0, y0):
    return ellipse(R, R, 0, x0, y0)


def build_incl(fem, r):
    points = circle(r, 0, 0)
    fem.make_inclusion(points)


def matprint(mat):
    print("")
    for x in mat:
        print("{:.3f}    {:.3f}  ".format(x[0], x[1]))
    print("")


def init_pattern():
    # # material pattern
    nmatx, nmaty = 2 ** 8, 2 ** 8
    mat = genmat.MaterialDensity()
    mat.sym8 = True
    mat.n_x, mat.n_y, mat.n_z = nmatx, nmaty, 1
    mat.p_seed = mat.mat_rand
    mat.nb_threshold = 2
    # mat.sym8 = True
    sigma_pat = 20
    mat.ratio_filter = [sigma_pat, sigma_pat, 1]
    return mat


def make_pattern(f, choice="rand", isample=0):
    mat = init_pattern()
    if choice is "circ":
        r = (f / pi) ** (1 / 2)
        pattern = genmat.ell_shapes(
            mat.mat_grid, rloc=[0, 0, 0], rwidth=[r, r, 0.1], m=2
        ).astype(float)
    elif choice is "one_dim":
        pattern = mat.discrete_pattern
        ps = pattern.shape
        pattern = np.zeros_like(pattern)
        pattern[0 : int(ps[0] * f), :, :] = 0
    elif choice is "rand":
        name = "random_samples_f_" + "{:.1f}".format(f)
        npzfile = np.load(os.path.join(mat_folder, name + ".npz"))
        sample = npzfile["sample"]
        pattern = sample[isample, :, :]
        pattern = pattern.reshape((pattern.shape[0], pattern.shape[1], 1))
    elif choice is "pat":
        pattern = mat.discrete_pattern
    else:
        raise TypeError
    mat.pattern = pattern
    return mat


def init_es(f, E_bias, incl=True, mat=None):
    r = (f / pi) ** (1 / 2)
    #####################################
    # # electrostatics
    #####################################
    fem_es = model_es.FemModel()

    fem_es.gmsh_verbose = 0  #: str: Gmsh verbose (int between 0 and 4)
    fem_es.getdp_verbose = 0  #: str: GetDP verbose (int between 0 and 4)
    fem_es.python_verbose = 0
    #: str: GetDP verbose (int between 0 and 1)
    fem_es.parmesh = 11
    fem_es.parmesh_incl = 11
    fem_es.type_des = "elements"
    fem_es.inclusion_flag = incl
    # if not incl:
    # fem_es.quad_mesh_flag=True
    fem_es.dx = 1
    fem_es.dy = 1
    fem_es.dz = 1
    fem_es.r = r
    fem_es.f = f
    fem_es.eps_incl = eps_incl
    fem_es.eps_host = eps_incl
    fem_es.E_static = E_bias
    fem_es.coupling_flag = True
    fem_es.tmp_dir = "./estat"
    fem_es.tmp_dir = tempfile.mkdtemp(prefix="/tmp/benjaminv.")
    fem_es.initialize()
    if incl:
        build_incl(fem_es, r)
        fem_es.make_mesh()
    else:
        if mat:
            fem_es.mat = mat
        else:
            fem_es.mat = make_pattern(f)

        fem_es.make_mesh()
        # fem_es.register_pattern(mat.pattern, mat.threshold_val)

    return fem_es


def init_hom(fem_es):
    #####################################
    # # homogenization
    #####################################
    fem_hom = model_hom.TwoScale2D()
    fem_hom.gmsh_verbose = fem_es.gmsh_verbose
    fem_hom.getdp_verbose = fem_es.getdp_verbose
    fem_hom.python_verbose = fem_es.python_verbose
    fem_hom.parmesh = fem_es.parmesh
    fem_hom.parmesh_incl = fem_es.parmesh_incl
    fem_hom.type_des = fem_es.type_des
    fem_hom.inclusion_flag = fem_es.inclusion_flag
    fem_hom.quad_mesh_flag = fem_es.quad_mesh_flag
    fem_hom.dx = fem_es.dx
    fem_hom.dy = fem_es.dy
    fem_hom.dz = fem_es.dz
    fem_hom.eps_incl = fem_es.eps_incl
    fem_hom.eps_host = fem_es.eps_host
    fem_hom.tmp_dir = tempfile.mkdtemp(prefix="/tmp/benjaminv.")
    fem_hom.coupling_flag = fem_es.coupling_flag
    fem_hom.aniso = True

    fem_hom.initialize()
    if fem_es.inclusion_flag:
        build_incl(fem_hom, fem_es.r)
        fem_hom.make_mesh()
    else:
        fem_hom.mat = fem_es.mat
        fem_hom.make_mesh()
        # fem_hom.register_pattern(mat.pattern, mat.threshold_val)
    return fem_hom


def compute_hom_pb(fem_hom, epsi, verbose=False):
    """Computes the homogenization problem.

    Parameters
    ----------
    fem_hom : instance of TwoScale2D class
        The two-scale FEM problem
    epsi : array-like
        The permittivity tensor diagonal elements, (xx, yy and zz components).
        Size is (3*nvar) where nvar is the number of points where the permittivity
        is defined in the FEM model (ie len(fem_hom.des[0])).

    Returns
    -------
    tuple (eps_hom, fem_hom)
        eps_hom: the permittivity tensor (3*3 numpy array)
        fem_hom: instance of TwoScale2D class

    """
    # interp = not fem_hom.inclusion_flag
    interp = False
    make_pos_tensor_eps(fem_hom, epsi, interp=interp)
    fem_hom.compute_solution()
    eps_hom = fem_hom.compute_epsilon_eff()
    if verbose:
        print("eps_hom")
        matprint(eps_hom)
    return eps_hom, fem_hom


def compute_elstat_pb(fem_es, epsi):
    # interp = not fem_es.inclusion_flag
    interp = False
    make_pos_tensor_eps(fem_es, epsi, interp=interp)
    fem_es.compute_solution()
    E = fem_es.postpro_electrostatic_field()
    E = np.real(E)
    return E, fem_es


def mat2des(fem):
    f = fem.make_fdens(fem.mat.pattern)
    Vdes = f(fem.des[1])
    return Vdes


def coupling_loop(fem_es, epsi, tol=1e-2, verbose=False):
    conv = False
    o = np.ones_like(epsi)
    Eold = fem_es.E_static * o, 0 * o, 0 * o
    while not conv:
        # if fem_es.inclusion_flag:
        E, fem_es = compute_elstat_pb(fem_es, epsi)
        epsi = epsilonr_ferroelectric(E.real, dc=True)
        if not fem_es.inclusion_flag:
            id = mat2des(fem_es)
            epsi_ = np.ones_like(E.real, dtype=complex)
            epsi_[:, id == 1] = eps_incl
            epsi_[:, id == 0] = epsi[:, id == 0]
            epsi = epsi_

        normdiff = normvec(E - Eold)
        normold = normvec(Eold)
        cv = np.mean(normdiff) / np.mean(normold)

        conv = cv < tol
        Eold = E
        if verbose:
            print("error: ", cv)
    # fem_es.postpro_fields(filetype="pos")
    # fem_es.open_gmsh_gui()
    return E, epsi, fem_es


def main(f, E_bias, coupling=True, incl=True, mat=None, verbose=True, rmtmpdir=True):
    fem_es = init_es(f, E_bias, incl=incl, mat=mat)
    nvar = len(fem_es.des[0])
    eps_f = epsilonr_ferroelectric(E_bias, dc=True)
    eps_f0 = epsilonr_ferroelectric(0, dc=True)
    id = np.ones(nvar)
    if incl:
        epsi = id * eps_f, id * eps_f0, id * eps_f0

    else:
        id = mat2des(fem_es)
        epsi_xx = np.ones_like(id, dtype=complex) * eps_incl
        epsi_yy = np.ones_like(id, dtype=complex) * eps_incl
        epsi_zz = np.ones_like(id, dtype=complex) * eps_incl
        epsi_xx[id == 0] = eps_f
        epsi_yy[id == 0] = eps_f0
        epsi_zz[id == 0] = eps_f0
        epsi = epsi_xx, epsi_yy, epsi_zz
    E = E_bias * id, 0 * id, 0 * id
    if E_bias != 0:
        if coupling:
            E, epsi, fem_es = coupling_loop(fem_es, epsi, verbose=verbose)
    fem_hom = init_hom(fem_es)
    epsi = epsilonr_ferroelectric(E)
    eps_hom, fem_hom = compute_hom_pb(fem_hom, epsi, verbose=verbose)
    if rmtmpdir:
        fem_hom.rm_tmp_dir()
        fem_es.rm_tmp_dir()
    return eps_hom, epsi, E, fem_hom, fem_es


def main_circle(params, save=False, coupling=True):
    E_bias, f = params
    print("Parameters: E = {:.2f}MV/m - f = {:.2f} ".format(E_bias, f))
    eps_hom, epsi, E, fem_hom, fem_es = main(f, E_bias, coupling=coupling)
    if save:
        fname = "circle_f_{:.2f}_E_{:.2f}".format(f, E_bias)
        if not coupling:
            fname += "_uncoupled"
        saveddir = os.path.join(data_folder, "circ_rods")
        try:
            os.mkdir(saveddir)
        except FileExistsError:
            pass
        np.savez(
            os.path.join(saveddir, fname + ".npz"),
            eps_hom=eps_hom,
            f=f,
            E_bias=E_bias,
            epsi=epsi,
            E=E,
        )


def main_circle_pattern(params, save=False, coupling=True):
    E_bias, f = params
    print("Parameters: E = {:.2f}MV/m - f = {:.2f} ".format(E_bias, f))
    mat = make_pattern(f, choice="circ")
    eps_hom, epsi, E, fem_hom, fem_es = main(
        f, E_bias, incl=False, coupling=coupling, mat=mat
    )
    if save:
        fname = "circle_f_{:.2f}_E_{:.2f}".format(f, E_bias)
        if not coupling:
            fname += "_uncoupled"
        saveddir = os.path.join(data_folder, "circ_rods_pattern")
        try:
            os.mkdir(saveddir)
        except FileExistsError:
            pass
        np.savez(
            os.path.join(saveddir, fname + ".npz"),
            eps_hom=eps_hom,
            f=f,
            E_bias=E_bias,
            epsi=epsi,
            E=E,
        )


def main_rand(params, save=False, coupling=True):
    E_bias, f = params
    print("Parameters: E = {:.2f}MV/m - f = {:.2f} ".format(E_bias, f))
    for isample in range(21):
        print("  sample {}".format(isample))
        mat = make_pattern(f, choice="rand", isample=isample)
        eps_hom, epsi, E, fem_hom, fem_es = main(
            f, E_bias, incl=False, coupling=coupling, mat=mat
        )
        if save:
            fname = "rand_f_{:.2f}_E_{:.2f}_sample_{}".format(f, E_bias, isample)
            if not coupling:
                fname += "_uncoupled"
            saveddir = os.path.join(data_folder, "rand_circ_rods")
            try:
                os.mkdir(saveddir)
            except FileExistsError:
                pass
            np.savez(
                os.path.join(saveddir, fname + ".npz"),
                eps_hom=eps_hom,
                f=f,
                E_bias=E_bias,
                epsi=epsi,
                E=E,
            )


def normvec(E):
    n = 0
    for e in E:
        n += np.abs(e) ** 2
    return np.sqrt(n)


f_parallel = parallel(main_circle, partype="scoop")
# f_parallel = parallel(main_circle_pattern, partype="scoop")

# f_parallel = parallel(main_rand, partype="scoop")

nE = 2
Ebias = np.linspace(0, 2, nE)
nF = 1
Fincl = np.linspace(0.1, 0.5, nF)
E1, F1 = np.meshgrid(Ebias, Fincl)
params = np.vstack((E1.ravel(), F1.ravel())).T


if __name__ == "__main__":

    f_parallel(params, save=True)
    # f_parallel(params, save=True, coupling=False)

    eps_f = epsilonr_ferroelectric(Ebias[-1], dc=False)
    eps_f0 = epsilonr_ferroelectric(0, dc=False)

    print(eps_f0 / eps_f)

#
#
# if __name__ == "__main__":
#     incl = False
#     f = 0.4
#     nbias = 7
#     bias = np.linspace(0, 2, nbias)
#     eps_f = epsilonr_ferroelectric(bias)
#     eps_f0 = epsilonr_ferroelectric(0)
#     eps_harm = harmonic_mean(f, eps_incl, eps_f)
#     eps_ar = arithmetic_mean(f, eps_incl, eps_f)
#
#     eps_c, eps_s = [], []
#     epsilon_c, epsilon_s = [], []
#     E_c, E_s = [], []
#
#     mat = make_pattern(f)
#     for E_bias in bias:
#         print("")
#         print("*" * 33)
#         print("Bias = {} MV/m".format(E_bias))
#         print("*" * 33)
#         print("")
#         print("Coupled --------------")
#         eps_hom, epsi, E, fem_hom, fem_es = main(
#             f, E_bias, coupling=True, incl=incl, mat=mat
#         )
#         print("eps_hom")
#         matprint(eps_hom)
#         print("Uncoupled ------------")
#         eps_hom_uc, epsi_uc, E_uc, fem_hom, fem_es = main(
#             f, E_bias, coupling=False, incl=incl, mat=mat
#         )
#         print("eps_hom_uc")
#         matprint(eps_hom_uc)
#         eps_c.append(eps_hom)
#         eps_s.append(eps_hom_uc)
#         epsilon_c.append(epsi)
#         epsilon_s.append(epsi_uc)
#         E_c.append(E)
#         E_s.append(E_uc)
#
#     eps_bias_c, eps_bias_s = (
#         np.zeros((2, nbias), dtype=complex),
#         np.zeros((2, nbias), dtype=complex),
#     )
#     eps_bias_c_mean, E_bias_c_mean = (
#         np.zeros((2, nbias), dtype=complex),
#         np.zeros((2, nbias), dtype=complex),
#     )
#     for i in range(2):
#         eps_bias_c[i, :] = np.array([_[i, i] for _ in eps_c])
#         eps_bias_s[i, :] = np.array([_[i, i] for _ in eps_s])
#         eps_bias_c_mean[i, :] = np.array([np.mean(_[i]) for _ in epsilon_c])
#         # E_bias_c_mean[i,:] = np.array([np.mean(_[i]) for _ in E_c])
#
#     eps_mg = [maxwell_garnett(f, eps_incl, eps_f)]
#     eps_mg.append(np.ones(nbias) * eps_mg[0][0])
#     eps_mg_mean = maxwell_garnett(f, eps_incl, eps_bias_c_mean)
#
#     # eps_ana_c, E_ana_c = [], []
#     # for E in bias:
#     #     eps_ = epsilonr_ferroelectric(E)
#     #     E, eps_ = coupling_loop_analytical(E, f, eps_incl, eps_)
#     #     eps_ana_c.append(eps_)
#     #     E_ana_c.append(E)
#     # eps_ana_c = np.array(eps_ana_c)
#     # eps_mg_ana_c = [maxwell_garnett(f, eps_incl, eps_ana_c)]
#     # eps_mg_ana_c.append(np.ones(nbias) * eps_mg_ana_c[0][0])
#
#     plt.close("all")
#     fig, ax = plt.subplots(2, figsize=(4, 8))
#     ax[0].plot(bias, eps_bias_c[0], "-r", label="coupled $\parallel$")
#     ax[0].plot(bias, eps_bias_s[0], "-b", label="uncoupled $\parallel$")
#     # ax[0].plot(bias, eps_mg[0], "-g", label="Maxwell-Garnett $\parallel$")
#     # ax[0].plot(bias, eps_mg_mean[0], "-m", label="numerical mean epsi MG $\parallel$")
#     # ax[0].plot(bias, eps_mg_ana_c[0], "-c", label="analytical mean epsi MG $\parallel$")
#
#     ax[0].legend()
#     ax[1].plot(bias, eps_bias_c[1], "-r", label="coupled $\perp$")
#     ax[1].plot(bias, eps_bias_s[1], "-b", label="uncoupled $\perp$")
#     # ax[1].plot(bias, eps_mg[1], "-g", label="Maxwell-Garnett $\perp$")
#     # ax[1].plot(bias, eps_mg_mean[1], "-m", label="numerical mean epsi MG $\perp$")
#     # ax[1].plot(bias, eps_mg_ana_c[1], "-c", label="analytical mean epsi MG $\perp$")
#
#     ax[1].legend()
#     ax[1].set_xlabel("Bias field $E$ (MV/m)")
#     ax[0].set_title("effective permittivity, $f = {}$".format(f))
#
#     ## loss tangent
#     tand_bias_c = losstan(eps_bias_c)
#     tand_bias_s = losstan(eps_bias_s)
#
#     fig, ax = plt.subplots(2, figsize=(4, 8))
#     ax[0].plot(bias, tand_bias_c[0], "-r", label="coupled $\parallel$")
#     ax[0].plot(bias, tand_bias_s[0], "-b", label="uncoupled $\parallel$")
#
#     ax[0].legend()
#     ax[1].plot(bias, tand_bias_c[1], "-r", label="coupled $\perp$")
#     ax[1].plot(bias, tand_bias_s[1], "-b", label="uncoupled $\perp$")
#
#     ax[1].legend()
#     ax[1].set_xlabel("Bias field $E$ (MV/m)")
#     ax[0].set_title("effective loss tangent, $f = {}$".format(f))
#
#     ## tunabilities
#
#     tunability_c = []
#     for e in eps_bias_c:
#         tunability_c.append(tunability(e, e[0]))
#     tunability_s = []
#     for e in eps_bias_s:
#         tunability_s.append(tunability(e, e[0]))
#
#     tunability_bulk = [tunability(eps_f, eps_f0)]
#     tbx = [1 for _ in range(nbias)]
#     tunability_bulk.append(tbx)
#     tunability_bulk.append(tbx)
#
#     fig, ax = plt.subplots(2, figsize=(4, 8))
#     ax[0].plot(bias, tunability_c[0], "-r", label="coupled $\parallel$")
#     ax[0].plot(bias, tunability_s[0], "-b", label="uncoupled $\parallel$")
#
#     ax[0].plot(bias, tunability_bulk[0], "-k", label="bulk $\parallel$")
#     ax[0].legend()
#     ax[1].plot(bias, tunability_c[1], "-r", label="coupled $\perp$")
#     ax[1].plot(bias, tunability_s[1], "-b", label="uncoupled $\perp$")
#
#     ax[1].plot(bias, tunability_bulk[1], "-k", label="bulk $\perp$")
#     ax[1].legend()
#     ax[1].set_xlabel("Bias field $E$ (MV/m)")
#     ax[0].set_title("tunability, $f = {}$".format(f))
#
#     # commutation quality factor
#
#     K_c = []
#     for e in eps_bias_c:
#         K_c.append(cqf(e, e[0]))
#     K_s = []
#     for e in eps_bias_s:
#         K_s.append(cqf(e, e[0]))
#
#     K_bulk = [cqf(eps_f, eps_f0)]
#     Ktmp = [0 for _ in range(nbias)]
#     K_bulk.append(Ktmp)
#     K_bulk.append(Ktmp)
#
#     fig, ax = plt.subplots(2, figsize=(4, 8))
#     ax[0].plot(bias, K_c[0], "-r", label="coupled $\parallel$")
#     ax[0].plot(bias, K_s[0], "-b", label="uncoupled $\parallel$")
#
#     ax[0].plot(bias, K_bulk[0], "-k", label="bulk $\parallel$")
#     ax[0].legend()
#     ax[1].plot(bias, K_c[1], "-r", label="coupled $\perp$")
#     ax[1].plot(bias, K_s[1], "-b", label="uncoupled $\perp$")
#
#     ax[1].plot(bias, K_bulk[1], "-k", label="bulk $\perp$")
#     ax[1].legend()
#     ax[1].set_xlabel("Bias field $E$ (MV/m)")
#     ax[0].set_title("commutation quality factor, $f = {}$".format(f))
