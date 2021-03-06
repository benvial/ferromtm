import numpy as np
from pytheas.material import genmat
from pytheas import periodic2D, utils
from pytheas.periodic2D import Periodic2D
from ferromtm.models.coupled2D import *
from ferromtm.visualization.plots import *
import tempfile

plt.close("all")


cv_dir_ = "circ_rods"
# cv_dir_ = "rand_circ_rods"
cv_dir = os.path.join(data_folder, cv_dir_, "convergence")


def load_arch(iter):
    fname = "cv_iter_{}.npz".format(iter)
    filename = os.path.join(cv_dir, fname)
    arch = np.load(filename)
    eps_hom = arch["eps_hom"]
    epsi = arch["epsi_map"]
    E = arch["E_map"]
    return epsi, E, eps_hom


pi = np.pi
f = 0.5
E0 = 2
rincl = (f / pi) ** 0.5
nincly = 5
dx = 1

ef = epsilonr_ferroelectric(E0)
epsi = [ef, 3]
matprop = [np.sqrt(e) for e in epsi]
fem = Periodic2D()

fem.lambda0 = 150  #: flt: incident wavelength

# opto-geometric parameters  -------------------------------------------
fem.d = 1  #: flt: period
fem.h_sup = fem.lambda0 * 2  #: flt: "thickness" superstrate
fem.h_sub = fem.lambda0 * 2  #: flt: "thickness" substrate
fem.h_layer1 = 1  #: flt: thickness layer 1
fem.h_layer2 = 1  #: flt: thickness layer 2
fem.h_des = nincly * dx  #: flt: thickness layer design
fem.h_pmltop = fem.lambda0  #: flt: thickness pml top
fem.h_pmlbot = fem.lambda0  #: flt: thickness pml bot
fem.a_pml = 1  #: flt: PMLs parameter, real part
fem.b_pml = 1  #: flt: PMLs parameter, imaginary part
fem.eps_sup = 1  #: flt: permittivity superstrate
fem.eps_sub = 1  #: flt: permittivity substrate
fem.eps_layer1 = 1  #: flt: permittivity layer 1
fem.eps_layer2 = 1  #: flt: permittivity layer 2
fem.eps_des = 1  #: flt: permittivity layer design

fem.theta_deg = 0.0  #: flt: incident angle
fem.pola = "TM"  #: str: polarization (TE or TM)
fem.lambda_mesh = 1  #: flt: incident wavelength
#: mesh parameters, correspond to a mesh size of lambda_mesh/(n*parmesh),
#: where n is the refractive index of the medium
fem.parmesh_des = 30
fem.parmesh = 0.22
fem.parmesh_pml = fem.parmesh * 2 / 3
fem.type_des = "elements"

debug = False
if debug:
    fem.gmsh_verbose = 4
    fem.getdp_verbose = 4
    fem.python_verbose = 1
else:
    fem.gmsh_verbose = 0
    fem.getdp_verbose = 0
    fem.python_verbose = 0

fem.aniso = True


def make_epsi(fem, epsi_map, nincly):
    epsi_xx = epsi_map[0]
    epsi_yy = epsi_map[1]
    n_x, n_y = epsi_xx.shape
    epsi_xx = np.tile(epsi_xx, (1, nincly))
    epsi_yy = np.tile(epsi_yy, (1, nincly))
    # plt.imshow(epsi_xx.real.T)

    epsi_xx = epsi_xx.reshape((n_x, n_y * nincly, 1))
    epsi_yy = epsi_xx.reshape((n_x, n_y * nincly, 1))
    id = np.ones_like(epsi_xx)

    # if homogenized:
    #     epsi_xx, epsi_yy = eps_hom[0, 0] * id, eps_hom[1, 1] * id

    epsi = epsi_xx, epsi_yy, id

    make_pos_tensor_eps(fem, epsi, interp=True)


def main_angle(t, coupled=True):
    if coupled:
        iter = 12  # periodic
    else:
        iter = 0
    epsi_map, E_map, eps_hom = load_arch(iter)
    fem.tmp_dir = tempfile.mkdtemp(prefix="/tmp/benjaminv.")
    fem.theta_deg = t
    ct = np.cos(t * pi / 180)
    fem.h_pmltop = fem.lambda0 / ct  #: flt: thickness pml top
    fem.h_pmlbot = fem.lambda0 / ct  #: flt: thickness pml bot
    fem.initialize()
    mesh = fem.make_mesh()
    make_epsi(fem, epsi_map, nincly)
    effs = []
    fem.compute_solution()
    effs = fem.diffraction_efficiencies()
    print("efficiencies", effs)
    return effs


def main_meta(t):
    return main_angle(t)


def main_meta_uncpl(t):
    return main_angle(t, coupled=False)


def rslab(h, lambda0, theta, epsxx, epsyy):
    """reflection coefficient for TM polarization of a
    slab with diagonal permittivity in vacuum"""

    t = pi * theta / 180
    ct = np.cos(t)
    st = np.sin(t)
    k0 = 2 * pi / lambda0
    k0x = k0 * st
    k0y = k0 * ct
    ky = k0 * (np.sqrt(epsxx * (1 - st ** 2 / epsyy)))

    r1 = (k0y - ky / epsxx) / (k0y + ky / epsxx)
    expo = np.exp(-2 * 1j * ky * h)
    r = r1 * (1 - expo) / (1 - r1 ** 2 * expo)
    return r


def main_hom(angle, coupled=True):
    if coupled:
        iter = 12  # periodic
    else:
        iter = 0
    _, _, eps_hom = load_arch(iter)
    r = rslab(fem.h_des, fem.lambda0, angle, eps_hom[0, 0], eps_hom[1, 1])
    return np.abs(r) ** 2


def get_key(out, k):
    return [e[k] for e in out]


def extract_effs(out):
    e = []
    for k in ["R", "T", "Q", "B"]:
        e.append(get_key(out, k))
    return e
