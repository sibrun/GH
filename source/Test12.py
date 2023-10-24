# fill nauty cache

import Parallel
import BufferedGeng


def fill_cache(ve):
    v, e = ve
    BufferedGeng.fill_cache(v, e, False)


if __name__ == "__main__":

    ve_pairs = [(v, l+v-1) for l in range(0, 11) for v in range(0,20)]
    # ve_pairs = [(v, l+v-1) for l in range(12, 13) for v in range(2*l-1)]
    Parallel.parallel(fill_cache, ve_pairs, n_jobs=4)
