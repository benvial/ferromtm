# -*- coding: utf-8 -*-
"""
Ploting a curve
==========================================

A dummy example
"""


import numpy as np
import matplotlib.pyplot as plt

##############################################################################
# Some text

x = np.linspace(0, 2 * np.pi, 300)
y = np.sin(x) ** 2

##############################################################################
# And finally plot it:

plt.close("all")
fig, ax = plt.subplots(1, figsize=(6, 4))
plt.plot(x, y, "r-", label=r"nice plot")
plt.xlabel(r"this is $x$")
plt.ylabel(r"this is $y$")
plt.legend(loc=0)
plt.show()
