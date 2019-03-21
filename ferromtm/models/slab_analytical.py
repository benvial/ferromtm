import numpy as np
from pytheas.material import genmat
from pytheas import periodic2D, utils
from pytheas.periodic2D import Periodic2D
from ferromtm.models.coupled2D import *
from ferromtm.visualization.plots import *
import tempfile

plt.close("all")


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


#
#
#
# def rslab_iso(h, lambda0, theta, eps):
#     """reflection coefficient for TM polarization of a
#     slab with diagonal permittivity in vacuum"""
#
#     t = pi * theta / 180
#     ct = np.cos(t)
#     st = np.sin(t)
#     k0 = 2 * pi / lambda0
#
#     k1 = k0
#     k2 = k0 * np.sqrt(eps)
#     a1 = -k1 * st
#
#     b1 = k1 * ct
#     b2 = np.sqrt(k2 ** 2 - a1 ** 2)
#
#     r12 = (b1 - b2 / eps) / (b1 + b2 / eps)
#
#     cti = np.conj(np.sqrt(1 - st ** 2 / eps + 0j))
#     deltai = cti * k0 * np.sqrt(eps) * h
#     expo = np.exp(-2 * 1j * deltai)
#
#     r = r12 * (1 - expo) / (1 - r12 ** 2 * expo)
#     t = 1 + r
#     return r, t


h = 5

lambda0 = 150

epsxx = 2.19112589e01 - 1.77562977e-01j
epsyy = 5.60140547e01 - 5.32664927e-01j

theta = np.linspace(0, 89, 111)


# lambda0 = np.linspace(0.1, 10, 100)

r = rslab(h, lambda0, theta, epsxx, epsyy)
R = np.abs(r) ** 2
plt.figure()
plt.plot(theta, R, "-r")
# plt.plot(lambda0, R)

epsxx = 2.98679870e01 - 2.72555535e-01j
epsyy = 5.40196497e01 - 5.11373051e-01j

ruc = rslab(h, lambda0, theta, epsxx, epsyy)

Ruc = np.abs(ruc) ** 2

plt.plot(theta, Ruc, "b-")
# plt.ylim((0, 1))
