import OrdinaryGraphComplex as OGC
import Profiling

reload(OGC)
reload(Profiling)

@Profiling.profile('log')
def ogc_main():
    skip_existing_files = False
    n_jobs = 1

    v_range = range(5, 10)
    l_range = range(5, 10)
    even_edges = True

    ogc = OGC.OrdinaryGC(v_range, l_range, even_edges)

    ogc.build_basis(skip_existing_files=skip_existing_files)
    ogc.build_operator_matrix(skip_existing_files=skip_existing_files, n_jobs=n_jobs)
    ogc.compute_ranks(skip_existing_files=skip_existing_files)
    ogc.compute_cohomology()
    ogc.store_cohomology_dim()
    ogc.plot_cohomology_dim()


if __name__ == "__main__":
    ogc_main()

    