import os
import subprocess
import numpy as np
import scipy as sc
from pytheas.tools import femio
from pytheas.basefem import BaseFEM, get_file_path

pi = np.pi


class FemModel(BaseFEM):
    """A class for computing electrostatic fields of
    a periodic 3D medium using a finite element model.
    See the base class :class:`BaseFEM` documentation for more info.
    """

    def __init__(
        self,
        #: flt: caracteristic length of the problem (typically the period)
        l_carac=1.0,
        #: flt: global mesh parameter
        #: `MeshElementSize = l_carac/(parmesh*n)`, `n`: refractive index
        parmesh=10.0,
        # opto-geometric parameters  -------------------------------------------
        dx=1,  #: flt: period x
        dy=1,  #: flt: period y
        dz=1,  #: flt: period z
        ax=0.25,  #: flt: ellipsoid principal axis length x
        ay=0.25,  #: flt: ellipsoid principal axis length y
        az=0.25,  #: flt: ellipsoid principal axis length z
        eps_host=1 - 0j,
        eps_incl=1 - 0j,
        dom_des=1002,  #: design domain number (check .geo/.pro files)
        dim=3,  #: dimension of the problem
        save_solution=False,
        type_des="nodes",
        inclusion_flag=False,
        coupling_flag=False,
        E_static=1,
    ):

        super().__init__()
        self.dir_path = get_file_path(__file__)

        #: flt: caracteristic length of the problem (typically the period)
        self.l_carac = l_carac

        #: flt: global mesh parameter
        #: `MeshElementSize = l_carac/(parmesh*n)`, `n`: refractive index
        self.parmesh = parmesh
        # opto-geometric parameters  -------------------------------------------
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.ax = ax
        self.ay = ay
        self.az = az
        self.eps_host = eps_host
        self.eps_incl = eps_incl
        self.dom_des = dom_des  #: design domain number (check .geo/.pro files)
        self.dim = dim  #: dimension of the problem

        self.save_solution = save_solution

        self.type_des = type_des
        self.lambda_mesh = l_carac
        self.eps_des = eps_host

        self.inclusion_flag = inclusion_flag
        self.E_static = E_static
        self.coupling_flag = coupling_flag

    celltype = "tetra"

    @property
    def corners_des(self):
        return -self.dx / 2, +self.dx / 2, -self.dy / 2, +self.dy / 2

    @property
    def domX_L(self):
        return -self.dx / 2

    @property
    def domX_R(self):
        return self.dx / 2

    @property
    def domY_B(self):
        return -self.dy / 2

    @property
    def domY_T(self):
        # + self.h_pmltop
        return self.dy / 2

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
        resolution = "electrostat"
        femio.solve_problem(
            resolution,
            self.path_pro,
            self.path_mesh,
            verbose=self.getdp_verbose,
            path_pos=self.path_pos,
            argstr=argstr,
        )

    def postpro_fields(self, filetype="txt"):
        self.print_progress("Postprocessing fields")
        self.postpro_choice("postop_fields", filetype)

    def postpro_electrostatic_field(self):
        self.print_progress("Postprocessing electrostatic field")
        subprocess.call(self.ppcmd("postop_E"))
        vect = femio.load_element_table_vect(self.tmp_dir + "/" + "Etot.txt")
        return np.array(vect)

    def get_field_map(self, name):
        field = femio.load_table(self.tmp_dir + "/" + name)
        return np.flipud(field.reshape((self.Niy, self.Nix))).T
