from ferromtm.models.slab import *

angle = list(np.linspace(0, 89, 100))
save_dir_ = "meta"
save_arch = os.path.join(data_folder, save_dir_, "efficiencies_angle.npz")


if __name__ == "__main__":
    from ferromtm.tools.parallelize import *

    main_meta_par = parallel(main_meta, partype="gridmap")
    main_meta_uncpl_par = parallel(main_meta_uncpl, partype="gridmap")
    out_meta = main_meta_par(angle)
    out_meta_uncpl = main_meta_uncpl_par(angle)
    np.savez(save_arch, out_meta=out_meta, out_meta_uncpl=out_meta_uncpl)
