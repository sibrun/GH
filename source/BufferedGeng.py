from sage.all import *
import OrdinaryGraphComplex
import OrdinaryMerkulovComplex
import OrdinaryGraphBiComplex
import HairyGraphComplex
import HairyGraphBiComplex
import BiColoredHairyGraphComplex
import BiColoredHairyGraphBiComplex
import Parameters
import StoreLoad
import LinboxInterface
import RheinfallInterface
import CHairyGraphComplex
import ForestedGraphComplex
import WRHairyGraphComplex
import os
import SymmetricGraphComplex
import NautyInterface

# this is a hack... should be geng, with symlink to bin"
geng_path = "/root/nauty27r3/gengL"


def fill_cache(n_vertices, n_edges, onlyonevi=True):
    if n_vertices <= 0 or n_edges <= 0 or 3 * n_vertices > 2 * n_edges or n_edges > n_vertices * (n_vertices - 1) / 2:
        return
    args = "Cd3" if onlyonevi else "cd3"

    filename = os.path.join(Parameters.geng_cachedir,
                            f"gengcache_{args}_{n_vertices}_{n_edges}.g6")
    if os.path.exists(filename):
        return

    StoreLoad.generate_path(filename)

    nauty_string = f"{geng_path} -{args} {n_vertices} {n_edges}:{n_edges} {filename}"
    print("Running nauty: ", nauty_string)
    NautyInterface.run_sys_cmd(nauty_string)


for v in range(8):
    for e in range(15):
        fill_cache(v, e, False)
