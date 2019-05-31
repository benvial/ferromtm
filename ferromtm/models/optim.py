from ferromtm.models.coupled2D import *
from pytheas.optim import *
from pytheas import femio
from pytheas.optim import topopt
from pytheas.optim import TopOpt
import tempfile
from aotomat.tools.plottools import *

# plt.close("all")

eps_interp = [1, 21]


# import autograd.numpy as np  # Thinly-wrapped version of Numpy
# from autograd import grad as grad_auto
# def taylor_sine(x):  # Taylor approximation to sine function
#     ca = np.mean(x)
#     return ca
#
# grad_sine = grad(taylor_sine)
# print( "Gradient of sin(pi) is", grad_sine(np.linspace(0,1,100)))


def init_pattern():
    # # material pattern
    nmatx, nmaty = 2 ** 8, 2 ** 8
    mat = genmat.MaterialDensity()
    mat.xsym = True
    mat.n_x, mat.n_y, mat.n_z = nmatx, nmaty, 1
    mat.p_seed = mat.mat_rand
    mat.nb_threshold = 2
    # mat.sym8 = True
    sigma_pat = 20
    mat.ratio_filter = [sigma_pat, sigma_pat, 1]
    return mat


mat = init_pattern()
mat.pattern = mat.normalized_pattern

parmesh = 21

fem_es = init_es(0, 0, incl=False, mat=mat, parmesh=parmesh, mesh_refine=False)


debug = False
if debug:
    fem_es.getdp_verbose = 4
    fem_es.gmsh_verbose = 4
    fem_es.python_verbose = 1


fem_es.quad_mesh_flag = True
# fem_es.gmsh_verbose = 4
fem_hom = init_hom(fem_es, tmp_dir="./tmp/test0")
fem_hom.pola = "TM"

# fem_hom.open_gmsh_gui()


# ##########################################
# #########  OPTIMIZATION PARAMETERS  ######
# ##########################################
to = TopOpt(fem_hom)
to.type_des = fem_hom.type_des
to.algorithm = topopt.nlopt.LD_MMA
# to.algorithm = topopt.nlopt.GN_DIRECT_L

to.typeopt = "min"  # type of optimization "min" or "max"
to.pmin = 0  # minimum value
to.pmax = 1  # maximum value
to.m = 1  # interpolation order eps=(eps_min-eps_max)*x^m-eps_min
to.ptol_rel = 1.0e-10
to.ftol_rel = 1.0e-5
to.stopval = 1e-12
to.maxeval = 20  # maximum of function evaluation
to.Nitmax = 8  # maximum number of global iterations
to.N0 = 0  # initial global iterations
# to.beta = 1  # projection parameter
lmin = 1
to.rfilt = 0.0 * lmin  # filter radius
to.filt_weight = "gaussian"
to.dg_dp = 0
to.eps_interp = eps_interp
to.log_opt = False
to.plotconv = False
to.force_xsym = False

to.dp = 1e-7
to.m = 1

ratio_hdes = 1
n_x = 37
n_y = 36  # int(n_x * ratio_hdes) +1
n_z = 1

to.n_x, to.n_y, to.n_z = n_x, n_y, n_z


def compute_hom_pb_y(fem_hom, epsi, verbose=False):
    interp = not fem_hom.inclusion_flag
    interp = False
    make_pos_tensor_eps(fem_hom, epsi, interp=interp)
    fem_hom.y_flag = True

    # fem_hom.path_pos += " ./tmp/test0/source_adj.pos"
    # print(fem_hom.path_pos)
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


def get_deq_deps(self, interp_method="nearest"):
    deq_deps = self.fem.get_deq_deps()
    x_grid, y_grid = self.grid

    deq_deps_x, deq_deps_y = deq_deps
    deq_deps_x = self.mesh2grid(deq_deps_x, interp_method=interp_method)
    deq_deps_x_x = np.gradient(deq_deps_x.T)[0] / np.gradient(x_grid)[0]

    deq_deps_y = self.mesh2grid(deq_deps_y, interp_method=interp_method)
    deq_deps_y_y = np.gradient(deq_deps_y.T)[1] / np.gradient(y_grid)[0]

    deq_deps = deq_deps_x_x.T + deq_deps_y_y.T

    # deq_deps = deq_deps_y_y.T

    deq_deps_re = self.grid2mesh(deq_deps.real)
    deq_deps_im = self.grid2mesh(deq_deps.imag)
    deq_deps = deq_deps_re + 1j * deq_deps_im
    return deq_deps


def get_sensitivity(to, p, filt=True, proj=True, interp_method="cubic"):
    to.fem.print_progress("Retrieving sensitivity")
    adjoint = to.get_adjoint()
    epsilon, depsilon_dp = to.make_epsilon(p, filt=filt, proj=proj, grad=True)
    deq_deps = get_deq_deps(to, interp_method=interp_method) * (-1 / epsilon ** 2)
    # deq_deps = to.get_deq_deps(interp_method=interp_method) * (-1 / epsilon ** 2)
    deq_deps = 1
    sens = to.dg_dp + np.real(adjoint * deq_deps * depsilon_dp)
    # plt_field(depsilon_dp, title="depsilon_dp")
    # plt_field(adjoint, title="adjoint")
    # plt_field(deq_deps, title="deq_deps")
    # plt_field(sens, title="sensitivities")
    return sens


def plt_field(a, title=None):
    plt.figure()
    plt.clf()
    a = to.mesh2grid(a.real, interp_method="nearest")
    plt.imshow(a)
    plt.colorbar()
    if title:
        plt.title(title)
    plt.pause(0.1)


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
    retgrad=False,
):
    sens_ana = np.size(grad) > 0
    fem_hom.adjoint = sens_ana
    epsilon, depsilon_dp = to.make_epsilon(p, filt=filt, proj=proj, grad=True)
    epsi = epsilon, epsilon, epsilon
    eps_hom_xx, fem_hom = compute_hom_pb_y(fem_hom, epsi, verbose=verbose)
    print("eps_hom_xx = ", eps_hom_xx)
    # obj0 = to.get_objective()
    eps_obj = 7
    obj = np.abs(1 / eps_obj - 1 / eps_hom_xx) ** 2 * eps_obj ** 2

    if sens_ana:
        sens = get_sensitivity(to, p, filt=filt, proj=proj)
    else:
        sens = 0
    # #
    # fem_hom.postpro_fields(filetype="pos")
    # fem_hom.open_gmsh_gui()

    # adj = to.get_adjoint()
    # plt_field(adj)

    if rmtmpdir:
        fem_hom.rm_tmp_dir()
        fem_es.rm_tmp_dir()
    print(("   objective =  %s " % obj))
    print("-" * 44)

    to.param_history.append(p)
    to.tot_obj_history.append(obj)
    to.Nit_tot += 1
    to.Nit_loc += 1

    if to.plotconv:
        make_plots(to, p, filt=filt, proj=proj)

    grad[:] = sens
    if retgrad:
        return obj, grad
    else:
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
        opt_f = f_obj(popt, np.array([]), filt=False, proj=False)
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
    np.random.seed(22)
    mat.p_seed = np.random.random(mat.pattern.shape)
    p0 = to.random_pattern(mat)
    print(len(p0))

    # out = main_opt(p0)

    fd = False

    if fd:
        grad_fd = []
        f = f_obj(p0, np.array([]))
        dp = 1e-3
        for i, p in enumerate(p0):
            print("iteration {0}/{1}".format(i, len(p0)))
            p1 = np.copy(p0)
            p1[i] += dp
            df = f_obj(p1, np.array([]))
            fd = (df - f) / dp
            print(fd)
            grad_fd.append(fd)
        grad_fd = np.array(grad_fd)
        np.savez("test_grad_fd.npz", grad_fd=grad_fd)
    else:
        grad_fd = np.load("test_grad_fd.npz")["grad_fd"]

    f, grad_adj = f_obj(p0, np.ones(len(p0)), retgrad=True)

    plt_field(grad_fd, "sens fd")

    # plt_field(grad_adj, "sens adj")

    adjoint = to.get_adjoint()
    epsilon, depsilon_dp = to.make_epsilon(p0, filt=True, proj=True, grad=True)
    deq_deps = to.get_deq_deps(interp_method="cubic") * (-1 / epsilon ** 2)
    # deq_deps = get_deq_deps(to, interp_method="cubic")* (epsilon )
    deq_deps = 1
    sens = to.dg_dp + np.real(adjoint * deq_deps * depsilon_dp)
    # plt_field(depsilon_dp, title="depsilon_dp")
    plt_field(adjoint, title="adjoint")
    # plt_field(deq_deps, title="deq_deps")
    plt_field(sens, title="sens adj")
    # plt_field(epsilon, title="epsilon")
    plt_field((-1 / epsilon ** 2), title="(-1/epsilon ** 2)")

    #
    # test =  adjoint * deq_deps
    # plt_field(test, title="sensitivities")

    # p0 = 0.5*np.ones_like(p0)
    # p0=np.random.random(len(p0))
    # Ebias = 0
    # c = f_obj(p0, f_obj,coupling=True, rmtmpdir=False)
    # out = main_opt(p0)
