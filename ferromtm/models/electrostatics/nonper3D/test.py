#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


import importlib
import femmodel

importlib.reload(femmodel)
from femmodel import *


from ferromtm.tools.utils import *
from ferromtm.models.theo import *


os.system("killall gmsh")

fem = Estat3D()
fem.rm_tmp_dir()


eps = 500

fem.debug = False

fem.eps_des = eps
fem.eps_electrode = 1  # -1100j
fem.eps_incl = 2.4
fem.parmesh_des = 5
fem.parmesh_host = 5
fem.parmesh_incl = 5
fem.parmesh_gap = 5
fem.parmesh_electrodes = 5

fem.lx_el = 6
fem.ly_el = 4.4
fem.lz_el = 0.1
fem.gap = 3.8


fem.hx_des = fem.lx_el * 2 + fem.gap + 3
fem.hy_des = fem.ly_el + 3
fem.hz_des = 0.6

fem.hx_box = fem.hx_des + 4
fem.hy_box = fem.hy_des + 4
fem.hz_box = fem.hz_des + fem.lz_el + 3


fem.R_hole = 1.1 / 2
fem.d_hole = 0.3 + fem.R_hole * 2


fem.Ebias = 2
fem.el_order = 1

fem.initialize()
#
# fem.open_gmsh_gui()
# cds

fem.make_mesh()
nvar = len(fem.des[0])

#
# fem.compute_solution()
# fem.postpro_fields_pos()

# os.system("cp ./base/geometry.geo.opt ./tmp/geometry.geo.opt")
# fem.open_gmsh_gui()
#
# E, P = fem.postpro_mean_fields()
# eps_eff = P / E
# print("effective permittivity = {}".format(eps_eff))
#
# eps_new = eps / 2
#
# fem.eps_des = eps_new
#
# fem.compute_solution()
# E, P = fem.postpro_mean_fields()
# eps_eff_new = P / E
# print("effective permittivity = {}".format(eps_eff_new))
#
#
# print("tunability theory = {}".format(eps / eps_new))
# print("tunability effective = {}".format(eps_eff / eps_eff_new))

fem.dom_des = 2
info_des = fem.get_mesh_info()
nvar_des = len(info_des[2][0])
fem.dom_des = 6
info_gap = fem.get_mesh_info()
nvar_gap = len(info_gap[2][0])

info = info_des, info_gap
name = "eps_des_", "eps_gap_"


def couple(Ebias):
    # fem.eps_des = epsilonr_ferroelectric(fem.Ebias)
    # dom_des

    i = 0
    fem.coupling_flag = True
    while True:
        if i == 0:
            Edes = (
                np.ones(nvar_des) * Ebias,
                np.ones(nvar_des) * 0,
                np.ones(nvar_des) * 0,
            )
            Egap = (
                np.ones(nvar_gap) * Ebias,
                np.ones(nvar_gap) * 0,
                np.ones(nvar_gap) * 0,
            )
            E = Edes, Egap
        for j in range(2):
            epsixx = epsilonr_ferroelectric(E[j][0].real)
            epsiyy = epsilonr_ferroelectric(E[j][1].real)
            epsizz = epsilonr_ferroelectric(E[j][2].real)
            epsi = epsixx, epsiyy, epsizz
            # e = np.ones_like(fem.des[0]) * epsilonr_ferroelectric(fem.Ebias)

            fem.des = info[j][2]
            make_pos_tensor_eps(fem, epsi, interp=False, basename=name[j])

        fem.compute_solution()

        # fem.postpro_fields_pos()
        # fem.open_gmsh_gui()
        # xsa
        Emean, Pmean = fem.postpro_mean_fields()

        eps_hom_xx = Pmean[0] / Emean[0]
        print("eps_hom_xx = ", eps_hom_xx)
        E = fem.postpro_electrostatic_field()

        if i == 0:
            eps_hom_xx_no_coupling = np.copy(eps_hom_xx)
        if i > 0:
            cv = np.abs(1 - eps_hom_xx / eps_hom_xx_)
            print("  cv = ", cv)
            if cv < 1e-2:
                break

        eps_hom_xx_ = np.copy(eps_hom_xx)
        i += 1
    return eps_hom_xx, eps_hom_xx_no_coupling


Ebias = np.linspace(2, 2, 1)
eps_hom_xx = []
eps_hom_xx_no_coupling = []
for fem.Ebias in Ebias:
    print("-------------------")
    eps = couple(fem.Ebias)
    eps_hom_xx.append(eps[0])
    eps_hom_xx_no_coupling.append(eps[1])


from aotomat.tools.plottools import *

plt.close("all")

col1 = "#6078cf"
col2 = "#c13f3f"

plt.figure()
# plt.plot(E_Theo, eps_Theo, "s", color=col1, alpha=0.3, label="bulk meas")
plt.plot(
    E_Theo_h, eps_Theo_h * eps_Theo_h_0, "o", color=col2, alpha=0.3, label="mtm meas"
)
# plt.plot(Ebias, epsilonr_ferroelectric(Ebias)/epsilonr_ferroelectric(Ebias)[0],label="bulk", color=col1)
plt.plot(Ebias, eps_hom_xx, label="coupling", color=col2)
plt.plot(Ebias, eps_hom_xx_no_coupling, "--", label="no coupling", color=col2)
plt.xlabel("$E$ (kV/mm)")
plt.ylabel("relative permittivity")
plt.legend()


fem.postpro_fields_pos()
fem.open_gmsh_gui()
