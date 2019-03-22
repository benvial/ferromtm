import tetrachotomy
from slab_modes import *

# import importlib
# importlib.reload(tetrachotomy)
tetrachotomy.plot_rect = True
tetrachotomy.plot_circ = True
tetrachotomy.plot_poles = True

plt.close("all")
pi = np.pi
tols = (1e-8 * (1 + 1j), 1e-8 * (1 + 1j), 1e-8 * (1 + 1j))
par_integ = (1e-8, 1e-8, 13)
tol_pol = 1e-8 * (1 + 1j)
tol_res = 1e-8 * (1 + 1j)
inv_golden_number = 2 / (1 + np.sqrt(5))
ratio = inv_golden_number
ratio_circ = 1 - inv_golden_number
nref_max = 100
ratio_re, ratio_im = ratio, ratio


z0 = 0.02 - 0.001 * 1j
z1 = 0.1 + 0.005 * 1j

z0 = -1 - 1 * 1j
z1 = 1 + 1 * 1j

x0, x1, y0, y1 = z0.real, z1.real, z0.imag, z1.imag
#####################################

E_bias = 2
epsi_map, eps_hom = load_arch(E_bias, coupling=True)


def main_hom(fnorm, coupled=True):
    fnorm_ = fnorm.real / 10 + 1j * fnorm.imag / 100
    lambda0 = fem.d / fnorm_
    r = rslab(fem.h_des, lambda0, fem.theta_deg, eps_hom[0, 0], eps_hom[1, 1])
    return np.abs(r) ** 2


def func(z):
    return 1 / main_hom(z, coupled=True)


# #
x_ = np.linspace(z0.real, z1.real, 201)
y_ = np.linspace(z0.imag, z1.imag, 200)
x, y = np.meshgrid(x_, y_)

z = x + 1j * y
u = func(z)

plt.figure()
uplot = np.log10(np.abs(func(z)))
plt.pcolor(x_, y_, uplot)
plt.colorbar()

fig = plt.figure()
ax = fig.add_subplot(111)  # , aspect='equal')
ax.set_xlim((x0, x1))
ax.set_ylim((y0, y1))
plt.xlabel(r"Re $z$")
plt.ylabel(r"Im $z$")

poles, residues, nb_cuts = tetrachotomy.pole_hunt(
    func,
    z0,
    z1,
    tols=tols,
    ratio_re=ratio_re,
    ratio_im=ratio_re,
    nref_max=nref_max,
    ratio_circ=ratio_circ,
    tol_pol=tol_pol,
    tol_res=tol_res,
    par_integ=par_integ,
    poles=[],
    residues=[],
    nb_cuts=0,
    verbose=1,
)
print("poles = ", poles)
print("residues = ", residues)
print("nb_cuts = ", nb_cuts)


fig = plt.figure()
fig = plt.gcf()
ax = fig.add_subplot(111)  # , aspect='equal')
ax.set_xlim((x0, x1))
ax.set_ylim((y0, y1))
plt.xlabel(r"Re $z$")
plt.ylabel(r"Im $z$")
plt.gca().plot(np.real(poles), np.imag(poles), "or")
