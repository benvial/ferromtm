#!/usr/bin/env python


from ferromtm.models.coupled2D import *
from ferromtm.models.electrostatics.per2D import femmodel
import importlib

importlib.reload(femmodel)

from ferromtm.tools.utils import *


fem = femmodel.FemModel()
fem.rm_tmp_dir()

# fem.gmsh_verbose=4
# fem.getdp_verbose=4

fem.parmesh = 4
fem.parmesh_des = 4

fem.parmesh_incl = 9


fem.E_static = 1

epsi0 = epsilonr_ferroelectric(0)
epsiE = epsilonr_ferroelectric(fem.E_static)

fem.eps_host = epsiE

# fem.eps_incl = epsiE

fem.b_pml = 0
fem.h_pml = 0.1

fem.space2pml_L = 0.1
fem.space2pml_R = 0.1
fem.space2pml_B = 0.1
fem.space2pml_T = 0.1


fem.inclusion_flag = True

dx, dy = 1, 1
nb_inclx, nb_incly = 1, 1


fem.hx_des = nb_inclx * dx
fem.hy_des = nb_incly * dy
nb_incl = nb_inclx * nb_incly
fem.nb_incl = nb_incl


r = 0.4

# Rx = (0.1 + np.random.random(nb_incl) * 0.3) * dx
# Ry = (0.1 + np.random.random(nb_incl) * 0.3) * dy
# rot_ = np.random.random(nb_incl) * 2 * pi
x00 = fem.hx_des / 2 - dx / 2
y00 = fem.hy_des / 2 - dy / 2
X0 = np.linspace(-x00, x00, nb_inclx)
Y0 = np.linspace(-y00, y00, nb_incly)
X0, Y0 = np.meshgrid(X0, Y0)
X0 = X0.ravel()
Y0 = Y0.ravel()
# Y0 = np.ones(nb_incl) * 0
# Y0 = (-0.5 + np.random.random(nb_incl) * 1) * dy
Rx = np.ones(nb_incl) * r * dx
Ry = np.ones(nb_incl) * r * dy
rot_ = np.linspace(-pi / 2, pi / 2, nb_incl)
# rot_ = np.ones(nb_incl) * 0


fem.initialize()

if fem.inclusion_flag:
    i = 0
    for Rinclx, Rincly, rot_incl, x0, y0 in zip(Rx, Ry, rot_, X0, Y0):
        points = ellipse(Rinclx, Rincly, rot_incl, x0, y0)
        # fem.inclusion_filename_ = "ellipse{0}.geo".format(i)
        # fem.make_inclusion(points, startpoint=1000 * (i + 1))
        fem.make_inclusion(points, startpoint=1000)
        i += 1


fem.make_mesh()

# fem.open_gmsh_gui()

nvar = len(fem.des[0])
epsixx = np.ones(nvar) * epsiE
epsiyy = np.ones(nvar) * epsi0
epsizz = np.ones(nvar) * epsi0


epsi = epsixx, epsiyy, epsizz

make_pos_tensor_eps(fem, epsi, interp=False)
fem.compute_solution()
fem.postpro_fields(filetype="pos")
fem.open_gmsh_gui()
