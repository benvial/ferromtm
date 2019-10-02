# -*- coding: utf-8 -*-
# Author: Benjamin Vial
# License: MIT


import numpy as np

# from ..tools import femio
from pytheas.basefem import *


import geom
import importlib

importlib.reload(geom)
from geom import make_geom


class Estat3D(BaseFEM):
    """Finite element model of a 3D electrostatic problem 
        using Gmsh_ and GetDP_.

        .. _Gmsh:
            http://gmsh.info/
        .. _GetDP:
            http://getdp.info/
    """

    def __init__(
        self,
        analysis="direct",
        hx_des=5,  #: flt: thickness x design
        hy_des=3,  #: flt: thickness y design
        hz_des=0.3,  #: flt: thickness z design
        hx_box=7,  #: flt: thickness x box
        hy_box=4,  #: flt: thickness y box
        hz_box=1,  #: flt: thickness z box
        eps_des=3 - 0 * 1j,  #: flt: permittivity design
        eps_host=1 - 0 * 1j,  #: flt: permittivity host
        eps_electrode=12,
        eps_incl=2,
        lcar=1,  #: flt: characteristic length
        parmesh_des=6,
        parmesh_host=6,
        parmesh_gap=6,
        parmesh_incl=6,
        parmesh_electrodes=6,
        gap=3.8,
        lx_el=3,
        ly_el=4.4,
        lz_el=0.6,
        d_hole=1.2,
        R_hole=2,
        Ebias=1,
        el_order=1,
        coupling_flag=False,
    ):
        super().__init__()
        self.dir_path = get_file_path(__file__)
        self.analysis = analysis

        # opto-geometric parameters  -------------------------------------------
        self.eps_des = eps_des
        self.eps_host = eps_host
        self.eps_electrode = eps_electrode
        self.eps_incl = eps_incl
        self.hx_des = hx_des
        self.hy_des = hy_des
        self.hz_des = hz_des
        self.hx_box = hx_box
        self.hy_box = hy_box
        self.hz_box = hz_box
        self.R_hole = R_hole
        self.d_hole = d_hole
        self.gap = gap
        self.lx_el = lx_el
        self.ly_el = ly_el
        self.lz_el = lz_el
        self.lcar = lcar
        self.Ebias = Ebias

        self.parmesh_des = parmesh_des
        self.parmesh_host = parmesh_host
        self.parmesh_incl = parmesh_incl
        self.parmesh_electrodes = parmesh_electrodes

        self.coupling_flag = coupling_flag

        self.el_order = el_order

        self.bg_mesh = False

        # 2  #: design domain number (check .geo/.pro files)
        self.dom_des = 6

        # postprocessing -------------------------------------------------
        #: int: number of x integration points
        #: for postprocessing diffraction efficiencies
        self.ninterv_integ = 60
        #: int: number of z slices points
        #: for postprocessing diffraction efficiencies
        self.nb_slice = 3
        #: flt:  such that `scan_dist  = min(h_sup, hsub)/scan_dist_ratio`
        self.scan_dist_ratio = 5

        self.dim = 3

        self.adjoint = False
        self.recomb_des = True
        # self.celltype="tetra"

    @property
    def celltype(self):
        return "tetra"

    @property
    def corners_des(self):
        return (
            -self.hx_des / 2,
            +self.hx_des / 2,
            -self.hy_des / 2,
            +self.hy_des / 2,
            -self.hz_des / 2,
            +self.hz_des / 2,
        )

    def initialize(self, *args, **kwargs):
        make_geom()
        return super().initialize(*args, **kwargs)

    # def make_param_dict(self):
    #     param_dict = super().make_param_dict()
    #     return param_dict

    # def mesh_model(self):
    #     super().mesh_model(mesh_format="msh4")

    def compute_solution(self, **kwargs):
        res_list = ["estat"]
        return super().compute_solution(res_list=res_list)

    def postpro_epsilon(self):
        self.print_progress("Postprocessing permittivity")
        self.postprocess("postop_epsilon" + " -order 2")

    def postpro_mean_fields(self):
        self.print_progress("Postprocessing mean fields")
        subprocess.call(self.ppcmd("postop_mean"))
        E = femio.load_table_vect(self.tmp_dir + "/" + "int_Efield.txt")
        P = femio.load_table_vect(self.tmp_dir + "/" + "int_polarization.txt")
        return np.array(E), np.array(P)

    def postpro_electrostatic_field(self):
        self.print_progress("Postprocessing electrostatic field")
        # subprocess.call(self.ppcmd("postop_E"))
        # self.postprocess("postop_E" + " -order 2")
        self.postprocess("postop_E")
        edes = np.array(femio.load_element_table_vect(self.tmp_dir + "/" + "Edes.txt"))
        egap = np.array(femio.load_element_table_vect(self.tmp_dir + "/" + "Egap.txt"))
        return edes, egap
