import HairyGraphBiComplex as HGBC
import multiprocessing as mp

if __name__ == "__main__":

    gc = HGBC.HairyBiGC(range(6, 12), -3, True, False)

    gc.build_basis(n_jobs=1, info_tracker=True)

    gc.build_matrix(n_jobs=1, info_tracker=True)