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

        #: flt: incident plane wave wavelength in free space
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
        #: int: number of x points for postprocessing field maps
        self.Nix = 100
        self.Niy = 100

        self.eps_host = 1 - 0j
        self.eps_incl = 1 - 0j

        # opto-geometric parameters  -------------------------------------------
        self.dx = 1  #: flt: period x
        self.dy = 1  #: flt: period y

        self.switch = False

    @property
    def dom_des(self):
        if self.switch:
            return 2000
        else:
            return 1000

    @dom_des.setter
    def dom_des(self, dom_des):
        pass

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
        self._print_progress("Computing solution: electrostatic problem")
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

    def postpro_fields(self, filetype="txt"):
        self._print_progress("Postprocessing fields")
        self._postpro_choice("postop_fields", filetype)

    def postpro_electrostatic_field(self):
        self._print_progress("Postprocessing electrostatic field")
        subprocess.call(self._ppcmd("postop_E"))
        vect = femio.load_element_table_vect(self.tmp_dir + "/" + "Etot.txt")
        return np.array(vect)

    def get_field_map(self, name):
        field = femio.load_table(self.tmp_dir + "/" + name)
        return np.flipud(field.reshape((self.Niy, self.Nix))).T


if __name__ == "__main__":
    print("This is the femmodel module")
