from ferromtm.models.coupled2D import *
from pytheas.optim import *
from pytheas import femio
from pytheas.optim import topopt
from pytheas.optim import TopOpt
import tempfile
from aotomat.tools.plottools import *
from aotomat.tools.parallelize import *

plt.close("all")

eps_interp = [1, 21]


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

parmesh = 4

fem_es = init_es(0, 0, incl=False, mat=mat, parmesh=parmesh, mesh_refine=False)

debug = False
if debug:
    fem_es.getdp_verbose = 4
    fem_es.gmsh_verbose = 4
    fem_es.python_verbose = 1


fem_es.quad_mesh_flag = True
# fem_es.gmsh_verbose = 4
fem_hom = init_hom(fem_es)
fem_hom.pola = "TM"

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
to.stopval = 1e-6
to.maxeval = 3  # maximum of function evaluation
to.Nitmax = 1  # maximum number of global iterations
to.N0 = 0  # initial global iterations
# to.beta = 1  # projection parameter
lmin = 1
to.rfilt = 0.0 * lmin  # filter radius
to.filt_weight = "gaussian"
to.dg_dp = 0
to.eps_interp = eps_interp
to.log_opt = False
to.plotconv = True
to.force_xsym = False
to.threshold_final = False

to.dp = 1e-7
to.m = 1

ratio_hdes = 1
n_x = 111
n_y = 113  # int(n_x * ratio_hdes) +1
n_z = 1

to.n_x, to.n_y, to.n_z = n_x, n_y, n_z


def compute_hom_pb_y(fem_hom, epsi, verbose=False):
    make_pos_tensor_eps(fem_hom, epsi, interp=False)
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


def get_sensitivity(to, p, filt=True, proj=True, interp_method="cubic"):
    print(to.fem.tmp_dir)
    to.fem.print_progress("Retrieving sensitivity")
    adjoint = to.get_adjoint()
    epsilon, depsilon_dp = to.make_epsilon(p, filt=filt, proj=proj, grad=True)
    deq_deps = to.get_deq_deps(interp_method=interp_method) * (-1 / epsilon ** 2)
    sens = to.dg_dp + np.real(adjoint * deq_deps * depsilon_dp)
    # plt_field(depsilon_dp, title="depsilon_dp")
    # plt_field(adjoint, title="adjoint")
    # plt_field(deq_deps, title="deq_deps")
    # plt_field(sens, title="sensitivities")
    return sens


def plt_field(a, title=None):
    plt.figure()
    plt.clf()
    a = to.mesh2grid(a.real)
    plt.imshow(a)
    plt.colorbar()
    if title:
        plt.title(title)
    plt.pause(0.1)


def objectivefunc(p, fem_hom=fem_hom, filt=True, proj=True, verbose=False, init=True):
    if init:
        fem_hom = init_hom(fem_es)
        # print(fem_hom.tmp_dir)
    epsilon, depsilon_dp = to.make_epsilon(p, filt=filt, proj=proj, grad=True)
    epsi = epsilon, epsilon, epsilon
    eps_hom_xx, fem_hom = compute_hom_pb_y(fem_hom, epsi, verbose=False)
    if verbose:
        print("eps_hom_xx = ", eps_hom_xx)
    # obj0 = to.get_objective()
    eps_obj = 7
    obj = np.abs(1 / eps_obj - 1 / eps_hom_xx) ** 2 * eps_obj ** 2
    return obj, fem_hom


def f_obj(
    p,
    grad,
    verbose=False,
    rmtmpdir=False,
    filt=True,
    proj=True,
    fem_hom=fem_hom,
    retgrad=False,
    fd=False,
):
    p.setflags(write=1)
    p[np.isnan(p)] = 0.5
    print(fem_hom.tmp_dir)

    sens_ana = np.size(grad) > 0
    if not fd:
        fem_hom.adjoint = sens_ana

    obj, fem_hom = objectivefunc(
        p, fem_hom=fem_hom, filt=filt, proj=proj, verbose=True, init=fd
    )
    to.fem = fem_hom

    if sens_ana:
        if fd:
            sens = get_grad_fd(p, para=True)
        else:
            sens = get_sensitivity(to, p, filt=filt, proj=proj)

        # sens[np.isnan(sens)] = 0
        sens[np.isnan(sens)] = 0
    else:
        sens = 0

    print("p = \n", p)
    print("sens = \n", sens)
    # #
    # fem_hom.postpro_fields(filetype="pos")
    fem_hom.open_gmsh_gui()
    #
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


def main_opt(p0, fd=True):

    if fd:

        def f(p, grad):
            return f_obj(p, grad, fd=True, rmtmpdir=True)

    else:
        f = f_obj

    # ##### MAIN OPTIMIZATION LOOP ############
    popt, opt_f, opt = to.main_loop_topopt(f, p0)
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


def get_grad_fd(p0, para=False):
    f, _ = objectivefunc(p0)
    dp = 1e-6
    nvar = len(p0)
    P = np.zeros((nvar, nvar))
    df = []
    for i, p in enumerate(p0):
        # print("iteration {0}/{1}".format(i,len(p0)))
        p1 = np.copy(p0)
        p1[i] += dp
        P[i, :] = np.array(p1)
    if para:

        res = objectivefuncpara(P)
        df = [_[0] for _ in res]
    else:
        for p in P:
            df.append(objectivefunc(p))
    grad_fd = (np.array(df) - f) / dp
    return grad_fd


if __name__ == "__main__":
    fd = True
    # partype = "gridmap"
    partype = "multiprocessing"
    objectivefuncpara = parallel(objectivefunc, partype=partype)
    # define initial density p0
    np.random.seed(22)
    mat.p_seed = np.random.random(mat.pattern.shape)
    p0 = to.random_pattern(mat)
    print(len(p0))

    out = main_opt(p0, fd=fd)

    dasda

    #
    if fd:
        grad_fd = get_grad_fd(p0)

        np.savez("test_grad_fd.npz", grad_fd=grad_fd)
    else:
        grad_fd = np.load("test_grad_fd.npz")["grad_fd"]
    plt_field(grad_fd)

    f, grad_adj = f_obj(p0, np.ones(len(p0)), retgrad=True, rmtmpdir=False)
    plt_field(grad_adj)

    to.fem.print_progress("Retrieving sensitivity")
    adjoint = to.get_adjoint()
    epsilon, depsilon_dp = to.make_epsilon(p0, filt=True, proj=True, grad=True)
    deq_deps = to.get_deq_deps(interp_method="cubic") * (-1 / epsilon ** 2)
    sens = to.dg_dp + np.real(adjoint * deq_deps * depsilon_dp)
    # plt_field(depsilon_dp, title="depsilon_dp")
    plt_field(adjoint, title="adjoint")
    plt_field(deq_deps, title="deq_deps")
    plt_field(sens, title="sensitivities")
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
