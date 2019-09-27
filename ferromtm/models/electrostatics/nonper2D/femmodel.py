import os
import subprocess
import numpy as np
import scipy as sc
from pytheas.tools import femio
from pytheas.basefem import BaseFEM, get_file_path

pi = np.pi


class FemModel(BaseFEM):
    """A class for computing electrostatic fields of
    a periodic 2D medium using a finite element model.
    See the base class :class:`BaseFEM` documentation for more info.
    """

    def __init__(self):
        super().__init__()
        self.dir_path = get_file_path(__file__)

        #: flt: caracteristic length (for mesh)
        self.lambda0 = 1.0
        #: flt: static applied electric field
        self.E_static = 1.0
        self.theta = 0

        #: flt: global mesh parameter
        #: `MeshElementSize = lambda0/(parmesh*n)`, `n`: refractive index
        self.parmesh = 10.0
        self.quad_mesh_flag = False
        self.coupling_flag = True
        self.type_des = "elements"

        self.nb_incl = 1
        #: flt: wavelength to use for meshing
        self.lambda_mesh = 1.0

        # opto-geometric parameters  -------------------------------------------
        self.h_pml = 1.0  #: flt: thickness pml
        self.hx_des = 2.0  #: flt: x - thickness scattering box (design)
        self.hy_des = 2.0  #: flt: y - thickness scattering box
        self.a_pml = 1  #: flt: PMLs parameter, real part
        self.b_pml = 1  #: flt: PMLs parameter, imaginary part
        self.eps_host = 1 - 0j  #: flt: permittivity host
        self.eps_des = 1 - 0j  #: flt: permittivity scattering box
        self.eps_incl = 1 - 0j  #: flt: permittivity inclusion
        self.dom_des = 5  #: design domain number (check .geo/.pro files)

        # postprocessing -------------------------------------------------
        #: coords of point for PostProcessing
        self.xpp, self.ypp = 0, 0
        self.space2pml_L, self.space2pml_R = 1.0, 1.0
        self.space2pml_T, self.space2pml_B = 1.0, 1.0

        #: int: number of x points for postprocessing field maps
        self.Nix = 100
        self.Niy = 100

    @property
    def corners_des(self):
        return -self.hx_des / 2, +self.hx_des / 2, -self.hy_des / 2, +self.hy_des / 2

    @property
    def domX_L(self):
        return -self.hx_des / 2 - self.space2pml_L

    @property
    def domX_R(self):
        return self.hx_des / 2 + self.space2pml_R

    @property
    def domY_B(self):
        return -self.hy_des / 2 - self.space2pml_B

    @property
    def domY_T(self):
        return self.hy_des / 2 + self.space2pml_T

    def make_param_dict(self):
        param_dict = super().make_param_dict()
        # param_dict["save_solution"] = int(self.save_solution)
        param_dict["coupling_flag"] = int(self.coupling_flag)
        return param_dict

    def compute_solution(self, **kwargs):
        if self.pattern:
            self.update_epsilon_value()
        self.update_params()
        self.print_progress("Computing solution: electrostatic problem")
        argstr = "-petsc_prealloc 1500 -ksp_type preonly \
                 -pc_type lu -pc_factor_mat_solver_type mumps"
        resolution = "electrostat_scalar"
        femio.solve_problem(
            resolution,
            self.path_pro,
            self.path_mesh,
            verbose=self.getdp_verbose,
            path_pos=self.path_pos,
            argstr=argstr,
        )

    # def postpro_fields(self, filetype="txt"):
    #     self.print_progress("Postprocessing fields")
    #     self.postpro_choice("postop_fields", filetype)

    def postpro_electrostatic_field(self):
        self.print_progress("Postprocessing electrostatic field")
        subprocess.call(self.ppcmd("postop_E"))
        vect = femio.load_element_table_vect(self.tmp_dir + "/" + "Etot.txt")
        return np.array(vect)

    def get_field_map(self, name):
        field = femio.load_table(self.tmp_dir + "/" + name)
        return np.flipud(field.reshape((self.Niy, self.Nix))).T

    def postpro_mean_fields(self):
        self.print_progress("Postprocessing mean fields")
        subprocess.call(self.ppcmd("postop_mean"))
        E = femio.load_table_vect(self.tmp_dir + "/" + "int_Efield.txt")
        P = femio.load_table_vect(self.tmp_dir + "/" + "int_polarization.txt")
        return np.array(E), np.array(P)
