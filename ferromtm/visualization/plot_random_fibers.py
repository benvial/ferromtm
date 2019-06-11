#!/usr/bin/env python

from ferromtm.visualization.plots import *
from ferromtm import rootdir
import numpy as np
import os


if __name__ == "__main__":
    frange = np.array([0.1, 0.2, 0.3, 0.4, 0.5])
    dirsave = os.path.join(rootdir, "data", "mat")
    plt.close("all")
    for f in frange:
        name = "random_samples_f_" + "{0}.npz".format(f)
        npzfile = np.load(os.path.join(dirsave, name))
        sample = npzfile["sample"]
        Npltx, Nplty = 3, 7
        fig, axarr = plt.subplots(
            Npltx, Nplty, figsize=(Nplty, Npltx), sharex=True, sharey=True
        )
        fig.subplots_adjust(hspace=0.1, wspace=0.1)
        fig.suptitle(r"Random material samples, $f = $ " + str(f))

        for i, ax in enumerate(axarr.ravel()):
            ax.axis("off")
            b = sample[i, :, :]
            ax.imshow(b, cmap="Greens_r")

        figname = "random_samples_f_" + "{0}percent.eps".format(int(f * 100))
        plt.savefig(
            os.path.join(rootdir, "data", "figures", figname), bbox_inches="tight"
        )
        plt.pause(0.1)
