#!/usr/bin/env python


from ferromtm.models.coupled2D import *

# from ferromtm.models.electrostatics.per2D import femmodel
from ferromtm.models.electrostatics.nonper2D import femmodel
import importlib

importlib.reload(femmodel)

from ferromtm.tools.utils import *

from ferromtm.models.theo import *

from pytheas.homogenization import *
from theo import *

if __name__ == "__main__":

    fem = femmodel.FemModel()
    fem.rm_tmp_dir()

    fem.dom_des = 10

    # fem.gmsh_verbose=4
    # fem.getdp_verbose=4

    fem.parmesh = 11
    fem.parmesh_des = 11
    fem.parmesh_incl = 11

    fem.E_static = 1

    epsi0 = epsilonr_ferroelectric(0)
    epsiE = epsilonr_ferroelectric(fem.E_static)

    fem.eps_host = 1

    fem.eps_incl = 2.4

    fem.b_pml = 0
    fem.h_pml = 1

    fem.space2pml_L = 1
    fem.space2pml_R = 1
    fem.space2pml_B = 1
    fem.space2pml_T = 1

    fem.inclusion_flag = True

    r = 0.5 * 1.1
    dholes = 0.3 + r * 2

    dx, dy = 2 * r + dholes, 2 * r + dholes
    nb_inclx, nb_incly = 2, 3

    fem.hx_des = 3.8
    fem.hy_des = 4.4
    nb_incl = nb_inclx * nb_incly
    fem.nb_incl = nb_incl

    lboxx = fem.hx_des + fem.space2pml_L + fem.space2pml_R + 2 * fem.h_pml
    lboxy = fem.hy_des + fem.space2pml_T + fem.space2pml_B + 2 * fem.h_pml

    # Rx = (0.1 + np.random.random(nb_incl) * 0.3) * dx
    # Ry = (0.1 + np.random.random(nb_incl) * 0.3) * dy
    # rot_ = np.random.random(nb_incl) * 2 * pi
    x00 = -dholes / 2
    y00 = -dholes
    X0 = np.linspace(-x00, x00, nb_inclx)
    Y0 = np.linspace(-y00, y00, nb_incly)
    X0, Y0 = np.meshgrid(X0, Y0)
    X0 = X0.ravel()
    Y0 = Y0.ravel()
    # Y0 = np.ones(nb_incl) * 0
    # Y0 = (-0.5 + np.random.random(nb_incl) * 1) * dy
    Rx = np.ones(nb_incl) * r
    Ry = np.ones(nb_incl) * r
    rot_ = np.linspace(-pi / 2, pi / 2, nb_incl)
    # rot_ = np.ones(nb_incl) * 0

    fem.initialize()

    if fem.inclusion_flag:
        i = 0
        for Rinclx, Rincly, rot_incl, x0, y0 in zip(Rx, Ry, rot_, X0, Y0):
            # x0 += (1 - 2 * np.random.rand()) * Rinclx / 3
            # y0 += (1 - 2 * np.random.rand()) * Rincly / 3
            points = ellipse(Rinclx, Rincly, rot_incl, x0, y0)
            fem.inclusion_filename_ = "ellipse{0}.geo".format(i)
            fem.make_inclusion(points, startpoint=1000 * (i + 1))
            # fem.make_inclusion(points, startpoint=1000 )
            i += 1

    fem.make_mesh()

    # fem.open_gmsh_gui()

    nvar = len(fem.des[0])
    # epsixx = np.ones(nvar)*epsiE
    # epsiyy = np.ones(nvar)*epsi0
    # epsizz = np.ones(nvar)*epsi0
    # epsi = epsixx, epsiyy, epsizz

    def couple(E):
        i = 0
        while True:
            epsixx = epsilonr_ferroelectric(E[0].real)
            epsiyy = epsilonr_ferroelectric(E[1].real)
            epsizz = epsilonr_ferroelectric(E[2].real)
            epsi = epsixx, epsiyy, epsizz
            # e = np.ones_like(epsixx) * fem.eps_incl
            # epsi = e, e, e

            make_pos_tensor_eps(fem, epsi, interp=False)
            fem.compute_solution()
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

    eps_hom_xx = []
    eps_hom_xx_no_coupling = []

    Ebias = np.linspace(1e-5, 2, 11)
    # Ebias = np.linspace(2, 2, 1)

    S = fem.hx_des * fem.hy_des

    f = nb_incl * pi * r ** 2 / (S)
    eps_host = epsilonr_ferroelectric(Ebias)
    # eps_host = eps_Theo_0
    epsmg = maxwell_garnett(f, fem.eps_incl, eps_Theo_0, dim=2)

    print("Maxwell-Garnett = ", epsmg)
    print("eps_Theo_h_0 = ", eps_Theo_h_0)

    run = True
    save = False
    if run:
        for fem.E_static in Ebias:
            E = np.ones(nvar) * fem.E_static, np.ones(nvar) * 0, np.ones(nvar) * 0
            print("-------------------")
            eps = couple(E)
            eps_hom_xx.append(eps[0])
            eps_hom_xx_no_coupling.append(eps[1])
        if save:
            np.savez(
                "test_theo.npz",
                Ebias=Ebias,
                eps_hom_xx=eps_hom_xx,
                eps_hom_xx_no_coupling=eps_hom_xx_no_coupling,
            )

    else:
        arch = np.load("test_theo.npz")
        Ebias = arch["Ebias"]
        eps_hom_xx = arch["eps_hom_xx"]
        eps_hom_xx_no_coupling = arch["eps_hom_xx_no_coupling"]

    # fem.postpro_fields(filetype="pos")
    # fem.open_gmsh_gui()

    from aotomat.tools.plottools import *

    plt.close("all")

    col1 = "#6078cf"
    col2 = "#c13f3f"
    plt.figure()
    plt.plot(E_Theo, eps_Theo, "s", color=col1, alpha=0.3, label="bulk meas")
    plt.plot(E_Theo_h, eps_Theo_h, "o", color=col2, alpha=0.3, label="mtm meas")
    plt.plot(
        Ebias,
        epsilonr_ferroelectric(Ebias) / epsilonr_ferroelectric(Ebias)[0],
        label="bulk",
        color=col1,
    )
    plt.plot(Ebias, eps_hom_xx / eps_hom_xx[0], label="coupling", color=col2)
    plt.plot(
        Ebias,
        eps_hom_xx_no_coupling / eps_hom_xx_no_coupling[0],
        "--",
        label="no coupling",
        color=col2,
    )
    plt.xlabel("$E$ (kV/mm)")
    plt.ylabel("normalized permittivity")
    plt.legend()

    plt.figure()
    # plt.plot(E_Theo, eps_Theo, "s", color=col1, alpha=0.3, label="bulk meas")
    plt.plot(
        E_Theo_h,
        eps_Theo_h * eps_Theo_h_0,
        "o",
        color=col2,
        alpha=0.3,
        label="mtm meas",
    )
    # plt.plot(Ebias, epsilonr_ferroelectric(Ebias)/epsilonr_ferroelectric(Ebias)[0],label="bulk", color=col1)
    plt.plot(Ebias, eps_hom_xx, label="coupling", color=col2)
    plt.plot(Ebias, eps_hom_xx_no_coupling, "--", label="no coupling", color=col2)
    plt.xlabel("$E$ (kV/mm)")
    plt.ylabel("relative permittivity")
    plt.legend()
