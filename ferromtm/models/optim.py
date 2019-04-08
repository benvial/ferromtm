from ferromtm.models.coupled2D import *
from pytheas.optim import *
from pytheas import femio
from pytheas.optim import topopt
from pytheas.optim import TopOpt
import tempfile
from aotomat.tools.plottools import *

eps_interp = [1, 13]

mat = init_pattern()
mat.pattern = mat.normalized_pattern
parmesh = 14
fem_es = init_es(0, 0, incl=False, mat=mat, parmesh=parmesh, mesh_refine=False)
fem_hom = init_hom(fem_es)


# ##########################################
# #########  OPTIMIZATION PARAMETERS  ######
# ##########################################
to = TopOpt(fem_hom)
to.type_des = fem_hom.type_des
to.algorithm = topopt.nlopt.LD_MMA
to.typeopt = "min"  # type of optimization "min" or "max"
to.pmin = 0  # minimum value
to.pmax = 1  # maximum value
to.m = 1  # interpolation order eps=(eps_min-eps_max)*x^m-eps_min
to.ptol_rel = 1.0e-6
to.ftol_rel = 1.0e-12
to.stopval = None
to.maxeval = 10  # maximum of function evaluation
to.Nitmax = 8  # maximum number of global iterations
to.N0 = 0  # initial global iterations
# to.beta = 1  # projection parameter
lmin = 1
to.rfilt = 0.05 * lmin  # filter radius
to.filt_weight = "gaussian"
to.dg_dp = 0
to.eps_interp = eps_interp
to.log_opt = False

to.plotconv = True
to.force_xsym = False

to.dp = 1e-7
to.m = 1

ratio_hdes = 1
n_x = 101
n_y = 100  # int(n_x * ratio_hdes) +1
n_z = 1

to.n_x, to.n_y, to.n_z = n_x, n_y, n_z


def compute_hom_pb_y(fem_hom, epsi, verbose=False):
    interp = not fem_hom.inclusion_flag
    interp = False
    make_pos_tensor_eps(fem_hom, epsi, interp=interp)
    fem_hom.y_flag = True
    fem_hom.compute_solution()
    fem_hom.postprocessing()
    V = fem_hom.get_vol()
    phi_yy = femio.load_table(fem_hom.tmppath("Phiyy.txt")) / V
    int_inveps_yy = femio.load_table(fem_hom.tmppath("I_inveps_yy.txt")) / V
    eps_hom_xx = 1 / (int_inveps_yy + phi_yy)
    if verbose:
        print("int_inveps_yy = ", int_inveps_yy)
        print("phi_yy = ", phi_yy)
    return eps_hom_xx, fem_hom


def f_obj(
    p,
    grad,
    coupling=True,
    mat=mat,
    record_cv=False,
    verbose=False,
    rmtmpdir=False,
    parmesh=parmesh,
    sens_ana=False,
    filt=True,
    proj=True,
    fem_es=fem_es,
    fem_hom=fem_hom,
):
    fem_hom.adjoint = sens_ana
    epsilon = to.make_epsilon(p, filt=filt, proj=proj)
    epsi = epsilon, epsilon, epsilon
    eps_hom_xx, fem_hom = compute_hom_pb_y(fem_hom, epsi, verbose=verbose)
    obj = to.get_objective()
    # print("objective: ", obj)
    if sens_ana:
        sens = to.get_sensitivity(p, filt=filt, proj=proj)
    else:
        sens = 0
    # fem_hom.postpro_fields(filetype="pos")
    # fem_hom.open_gmsh_gui()
    if rmtmpdir:
        fem_hom.rm_tmp_dir()
        fem_es.rm_tmp_dir()
    print(("   objective =  %s " % obj))
    print("-" * 44)

    # to.obj_history.append(obj)
    # print(to.obj_history)
    to.param_history.append(p)
    to.tot_obj_history.append(obj)
    to.Nit_tot += 1
    to.Nit_loc += 1

    if to.plotconv:
        make_plots(to, p, filt=filt, proj=proj)
    return obj


def main_opt(p0):

    # ##### MAIN OPTIMIZATION LOOP ############
    popt, opt_f, opt = to.main_loop_topopt(f_obj, p0)
    print("optimum at ", popt)
    print("with value  = ", opt_f)
    print(popt)

    if to.threshold_final:
        print("\n")
        print("Final design")
        print("#" * 60)
        if to.force_xsym:
            popt = to.make_xsym(popt)
        popt_filt, _ = to.filter_param(popt, grad=False)
        popt = to.get_threshold_design(popt_filt)
        opt_f = f(popt, np.array([]), filt=False, proj=False)
        print("optimum at ", popt)
        print("with value  = ", opt_f)

    return popt, opt_f, to, fem_hom


def make_plots(to, p, filt=True, proj=True):
    # print("Plotting")
    epsilon = to.make_epsilon(p, filt=filt, proj=proj)
    qtplot = epsilon.real
    title = r"permittivity"
    to.plot_while_solving(
        qtplot,
        title=title,
        cmap="viridis",
        typeplot="interp",
        extent=(0, 1, 0, 1),
        interp_method="nearest",
    )


if __name__ == "__main__":
    # define initial density p0
    p0 = to.random_pattern(mat)
    Ebias = 0
    # c = f_obj(p0, f_obj,coupling=True, rmtmpdir=False)
    out = main_opt(p0)
