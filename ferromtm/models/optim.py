from ferromtm.models.coupled2D import *
from pytheas.optim import *
from pytheas import femio
from pytheas.optim import topopt
from pytheas.optim import TopOpt
import tempfile
from aotomat.tools.plottools import *

eps_interp = [1, 21]


# import autograd.numpy as np  # Thinly-wrapped version of Numpy
# from autograd import grad as grad_auto

#
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

parmesh = 11

fem_es = init_es(0, 0, incl=False, mat=mat, parmesh=parmesh, mesh_refine=False)

fem_es.quad_mesh_flag = False
# fem_es.gmsh_verbose = 4
fem_hom = init_hom(fem_es, tmp_dir="./tmp/test0")
fem_hom.pola = "TE"

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
n_x = 111
n_y = 111  # int(n_x * ratio_hdes) +1
n_z = 1

to.n_x, to.n_y, to.n_z = n_x, n_y, n_z


fem_hom.getdp_verbose = 0


def compute_hom_pb_y(fem_hom, epsi, verbose=False):
    interp = not fem_hom.inclusion_flag
    interp = False
    make_pos_tensor_eps(fem_hom, epsi, interp=interp)
    fem_hom.y_flag = True

    fem_hom.path_pos += " ./tmp/test0/source_adj.pos"
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
    sens_ana = np.size(grad) > 0
    fem_hom.adjoint = sens_ana

    epsilon, depsilon_dp = to.make_epsilon(p, filt=filt, proj=proj, grad=True)
    epsi = epsilon, epsilon, epsilon

    eps_hom_xx, fem_hom = compute_hom_pb_y(fem_hom, epsi, verbose=verbose)
    print("eps_hom_xx = ", eps_hom_xx)
    obj0 = to.get_objective()
    eps_obj = 18
    obj = np.abs(1 / eps_obj - obj0) ** 2 * (eps_obj) ** 2

    obj = obj0
    V = 1
    int_inveps_yy = femio.load_table(fem_hom.tmppath("I_inveps_yy.txt")) / V

    if sens_ana:
        sens = to.get_sensitivity(p, filt=filt, proj=proj) * int_inveps_yy.real
    else:
        sens = 0
    # fem_hom.postpro_fields(filetype="pos")
    # fem_hom.open_gmsh_gui()
    # ncdc

    # plt.clf()
    # sol = fem_hom.get_solution().real
    # # # print(adj)
    # # solplt = to.mesh2grid(sol.real)
    # # ux, uy = np.gradient(solplt)
    # # plt.imshow(ux)
    # # plt.colorbar()
    # # plt.pause(1)
    # #
    # x, y = to.grid
    # dx = x[1] - x[0]
    # dy = y[1] - y[0]
    # xsi = 1 / epsilon
    # xsi_tmp = to.mesh2grid(xsi)
    # u = to.mesh2grid(sol)
    #
    #
    def objective_func(u):
        ux, uy = np.gradient(u.T, edge_order=2)
        integ = xsi_tmp.T * (1 + uy / dy)
        xsihom = np.trapz(np.trapz(integ, y), x)
        # xsihom = np.mean(u)
        return xsihom

    #
    # def grad_objective_func(u, du=1e-2):
    #     g = np.zeros_like(u)
    #     N,M = u.shape
    #     f = objective_func(u)
    #     for ix in range(N):
    #         for iy in range(M):
    #             u_ = np.copy(u)
    #             u_[ix,iy] += du
    #             df = objective_func(u_)
    #             g[ix,iy]= (df-f)/du
    #     return g
    #
    def grad_objective_func(umesh):
        du = 1e-7
        df = []
        ugrid = to.mesh2grid(umesh)
        f = objective_func(ugrid)
        for i, u_ in enumerate(umesh):
            u_tmp = np.copy(umesh)
            u_tmp[i] += du
            ugrid = to.mesh2grid(u_tmp)
            df.append(objective_func(ugrid))
        df = np.array(df)
        # du = np.gradient(umesh)
        return (df - f) / du

    #

    # source_adj = grad_objective_func(sol)
    # fem_hom.make_eps_pos( fem_hom.des[0], -source_adj, posname="source_adj")

    # xsihom = objective_func(u)
    # print("eps_hom_xx test = ", 1 / xsihom)

    # grad_objective_func = grad_auto(objective_func)

    # dgdsol = grad_objective_func(sol)
    # dgdsol = to.mesh2grid(dgdsol)
    # plt.clf()
    # plt.imshow(dgdsol.real)
    # plt.colorbar()
    # plt.pause(2)
    # #
    # # dgdsol = grad_objective_func(u)
    # # plt.clf()
    # # plt.imshow(dgdsol.real)
    # # plt.colorbar()
    #

    # plt.clf()
    # # print(sens)
    # sensplt = to.mesh2grid(sens)
    # plt.imshow(sensplt)
    # plt.colorbar()
    # plt.pause(3)
    # #
    # plt.clf()
    # adj = to.get_adjoint()
    # # print(adj)
    # adjplt = to.mesh2grid(adj.real)
    # plt.imshow(adjplt)
    # plt.colorbar()
    # plt.pause(3)

    # plt.clf()
    # deq_deps = to.get_deq_deps()
    # # print(deq_deps)
    # deq_deps_plt = to.mesh2grid(deq_deps.real)
    # plt.imshow(deq_deps_plt)
    # plt.colorbar()
    # # # cds
    # #
    # plt.pause(3)

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

    grad[:] = sens
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
    # p0 = 0.5*np.ones_like(p0)
    # p0=np.random.random(len(p0))
    Ebias = 0
    # c = f_obj(p0, f_obj,coupling=True, rmtmpdir=False)
    out = main_opt(p0)
