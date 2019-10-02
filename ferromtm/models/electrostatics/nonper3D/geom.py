#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


import pygmsh
import os

# ext = False
def make_geom():
    geom_oc = pygmsh.opencascade.Geometry()

    geom_oc.add_raw_code('Include "parameters.dat";')
    geom_oc.add_raw_code("lc_des = lcar/(parmesh_des);")
    geom_oc.add_raw_code("lc_host = lcar/(parmesh_host);")
    geom_oc.add_raw_code("lc_holes = lcar/(parmesh_incl);")
    geom_oc.add_raw_code("lc_electrodes = lcar/(parmesh_electrodes);")
    geom_oc.add_raw_code("lc_gap = lcar/(parmesh_gap);")

    def make_des():
        des = geom_oc.add_rectangle(
            ["-hx_des/2", "-hy_des/2", "-hz_des/2"],
            "hx_des",
            "hy_des",
            char_length="lc_des",
        )
        _, des_ext, _ = geom_oc.extrude(des, [0, 0, "hz_des"])
        return des_ext

    def make_holes():
        rholes = ["R_hole"] * 6
        xholes = [
            "-d_hole / 2",
            "d_hole / 2",
            "-d_hole / 2",
            "d_hole / 2",
            "-d_hole / 2",
            " d_hole / 2",
        ]
        yholes = ["-d_hole", "0", "d_hole", "-d_hole", "0", "d_hole"]
        holes = []
        holes_ext = []
        for xh, yh, R_hole in zip(xholes, yholes, rholes):
            hole = geom_oc.add_disk([xh, yh, "-hz_des/2"], R_hole, char_length="lc_des")
            holes.append(hole)
            _, holes_ext_, _ = geom_oc.extrude(hole, [0, 0, "hz_des"])
            holes_ext.append(holes_ext_)
        return holes_ext

    def make_box():
        box = geom_oc.add_rectangle(
            ["-hx_box/2", "-hy_box/2", "-hz_box/2"],
            "hx_box",
            "hy_box",
            char_length="lc_host",
        )
        _, box, _ = geom_oc.extrude(box, [0, 0, "hz_box"])
        return box

    def make_gap():
        gap = geom_oc.add_rectangle(
            ["-gap/2", "-ly_el/2", "-hz_des/2"], "gap", "ly_el", char_length="lc_gap"
        )
        _, gap, _ = geom_oc.extrude(gap, [0, 0, "hz_des"])
        return gap

    def make_electrodes():
        electrodes = [None, None]
        electrode_left = geom_oc.add_rectangle(
            ["-gap / 2 - lx_el", "-ly_el / 2", "hz_des/2"],
            "lx_el",
            "ly_el",
            char_length="lc_electrodes",
        )
        _, electrodes[0], _ = geom_oc.extrude(electrode_left, [0, 0, "lz_el"])
        electrode_right = geom_oc.add_rectangle(
            ["gap / 2", "-ly_el / 2", "hz_des/2"],
            "lx_el",
            "ly_el",
            char_length="lc_electrodes",
        )
        _, electrodes[1], _ = geom_oc.extrude(electrode_right, [0, 0, "lz_el"])
        return electrodes

    # hole = geom_oc.add_disk(["0", "0", "-hz_des/2"], R_hole, char_length="lc_des")

    box = make_box()
    des = make_des()
    electrodes = make_electrodes()
    bkg = geom_oc.boolean_difference([box], [des] + electrodes)
    ###

    des = make_des()
    gap = make_gap()
    des = geom_oc.boolean_difference([des], [gap])

    gap = make_gap()
    holes = make_holes()
    gap = geom_oc.boolean_difference([gap], holes)

    ###
    electrodes = make_electrodes()
    holes = make_holes()

    printpoint = geom_oc.add_point([0, 0, 0])

    geom_oc.add_physical(bkg)
    geom_oc.add_physical(des)

    geom_oc.add_physical(holes)
    geom_oc.add_physical(electrodes[0])
    geom_oc.add_physical(electrodes[1])
    geom_oc.add_physical(gap)
    # geom_oc.add_physical([gap, des])
    geom_oc.add_physical(printpoint)

    geom_oc.add_raw_code("Coherence;")
    geom_oc.add_raw_code("Coherence;")
    geom_oc.add_raw_code("Coherence;")

    # geom_oc.add_raw_code("Mesh.Algorithm=6;")
    # geom_oc.add_raw_code("Mesh.Algorithm3D=4;")

    code = geom_oc.get_code().replace("'", "")

    dir_path = os.path.dirname(os.path.abspath(__file__))

    geo_filename = dir_path + "/base/geometry.geo"
    with open(geo_filename, "w") as f:
        f.write(code)


# import os
# dir_path = os.path.dirname(os.path.abspath(__file__))
#
# geo_filename = dir_path + "/base/geometry.geo"
# make_geom(ext=False)


# os.system("gmsh " + geo_filename)
