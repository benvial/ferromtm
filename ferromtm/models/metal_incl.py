#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


from ferromtm.models.coupled2D import *

E_bias = 2
f = 0.4
parmesh = 5
# eps_incl = 1000 - 1000j
# eps_incl = 3 - 0j


epsi_0 = epsilonr_ferroelectric(0)
epsi_E = epsilonr_ferroelectric(E_bias)

t = epsi_0.real / epsi_E.real

# EPS_INCL=[3]
EPS_INCL = np.linspace(1, 50, 15)


TC = []
TU = []

for eps_incl in EPS_INCL:

    parmesh = 5 / np.sqrt(eps_incl.real)

    eps_hom0, epsi0, E0, fem_hom0, fem_es0 = main(
        f=f,
        E_bias=0,
        coupling=False,
        eps_incl=eps_incl,
        parmesh=parmesh,
        rmtmpdir=False,
    )
    eps_hom, epsi, E, fem_hom, fem_es = main(
        f=f,
        E_bias=E_bias,
        coupling=False,
        eps_incl=eps_incl,
        parmesh=parmesh,
        rmtmpdir=False,
    )

    eps_hom_coupled, epsi_coupled, E_coupled, fem_hom_coupled, fem_es_coupled = main(
        f=f,
        E_bias=E_bias,
        coupling=True,
        eps_incl=eps_incl,
        parmesh=parmesh,
        rmtmpdir=False,
    )

    # fem_hom.postpro_fields_pos()

    # fem_hom.open_gmsh_gui()

    te = eps_hom0[0, 0].real / eps_hom[0, 0].real

    te_coupled = eps_hom0[0, 0].real / eps_hom_coupled[0, 0].real
    print("*" * 33)
    print("Tunability")
    print(f"  bulk {t}\n  uncoupled {te}\n  coupled {te_coupled}")
    print("*" * 33)

    TU.append(te)
    TC.append(te_coupled)


TB = np.ones_like(EPS_INCL) * t

from aotomat.tools.plottools import *

plt.plot(EPS_INCL, TB, "--k", label="bulk")
plt.plot(EPS_INCL, TU, c=aotomat_green, label="uncoupled")
plt.plot(EPS_INCL, TC, c=aotomat_orange, label="coupled")
plt.xlabel(r"inclusion permittivity $\varepsilon_{\rm incl}$")
plt.ylabel(r"tunability $n$")
plt.xlim((EPS_INCL[0], EPS_INCL[-1]))
plt.legend()
plt.tight_layout()
