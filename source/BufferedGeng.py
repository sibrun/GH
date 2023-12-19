
from sage.all import *
import Parameters
import StoreLoad
import os
import NautyInterface


# this is a hack... should be geng, with symlink to bin"
# geng_path = "/root/nauty27r3/gengL"
geng_path = "gengL"


def _get_geng_args_and_file(n_vertices, n_edges, onlyonevi=True):
    args = "Cd3" if onlyonevi else "cd3"
    filename = os.path.join(Parameters.geng_cachedir,
                            f"gengcache_{args}_{n_vertices}_{n_edges}.g6")
    return (args, filename)


def fill_cache(n_vertices, n_edges, onlyonevi=True):
    # WARNING: Ctrl+C-ing this computation results in an incomplete, i.e., corrupt, cache file!!!!!
    # Could be improved by using temp file then copying...
    if n_vertices <= 0 or n_edges <= 0 or 3 * n_vertices > 2 * n_edges or n_edges > n_vertices * (n_vertices - 1) / 2:
        return
    args, filename = _get_geng_args_and_file(n_vertices, n_edges, onlyonevi)
    if os.path.exists(filename):
        return

    StoreLoad.generate_path(filename)

    nauty_string = f"{geng_path} -{args} {n_vertices} {n_edges}:{n_edges} {filename}"
    print("Running nauty: ", nauty_string)
    NautyInterface.run_sys_cmd(nauty_string)


def list_simple_graphs_buffered(n_vertices, n_edges, onlyonevi=True):
    if n_vertices <= 0 or n_edges <= 0 or 3 * n_vertices > 2 * n_edges or n_edges > n_vertices * (n_vertices - 1) / 2:
        return []
    _, filename = _get_geng_args_and_file(n_vertices, n_edges, onlyonevi)

    if not os.path.exists(filename):
        print(f"Buffered geng: requested cache file {filename} does not exist. Generate cache first")
        raise ValueError(
            f"Buffered geng: requested cache file {filename} does not exist. Generate cache first")
    with open(filename, "r") as f:
        txt = f.read()
        if type(txt) is not str:
            txt = txt.decode("ascii")
        list_g6 = txt.splitlines()
        print(f"Buffered nauty generated {len(list_g6)} graphs.")

    return (Graph(g6) for g6 in list_g6)


# for l in range(12):
#     print(f"****** Loop order {l} ******")
#     for v in range(2*l-1):
#         e = l+v-1
#         fill_cache(v, e, False)
