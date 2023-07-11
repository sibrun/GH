"""This script generates the latex tables for the paper."""

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


# ***** only use if external hd
# Parameters.data_home_dir = "/Volumes/backup2/gh_data/"
# Parameters.data_dir = Parameters.data_home_dir + "data"
# Parameters.plots_dir = Parameters.data_home_dir + "plots"
# Parameters.ref_data_dir = Parameters.data_home_dir + "data_ref"
# Parameters.log_dir = Parameters.data_home_dir + "log"
##########

latexdir = os.path.join(Parameters.plots_dir, "latex", "paper")

latexfile_wrhairy_vs = os.path.join(latexdir, "wrhairy_vs.tex")
latexfile_wrhairy_ops = os.path.join(latexdir, "wrhairy_ops.tex")
latexfile_wrhairy_cohom = os.path.join(latexdir, "wrhairy_cohom.tex")

latexfile_ordinary_vs = os.path.join(latexdir, "ordinary_vs.tex")
latexfile_ordinary_ops = os.path.join(latexdir, "ordinary_ops.tex")
latexfile_ordinary_cohom = os.path.join(latexdir, "ordinary_cohom.tex")

latexfile_ordinaryme_vs = os.path.join(latexdir, "ordinaryme_vs.tex")
latexfile_ordinaryme_ops = os.path.join(latexdir, "ordinaryme_ops.tex")
latexfile_ordinaryme_cohom = os.path.join(latexdir, "ordinaryme_cohom.tex")

latexfile_hairy_vs = os.path.join(latexdir, "hairy_vs.tex")
latexfile_hairy_ops = os.path.join(latexdir, "hairy_ops.tex")
latexfile_hairy_cohom = os.path.join(latexdir, "hairy_cohom.tex")
latexfile_hairy_cohom_ee = os.path.join(latexdir, "hairy_cohom_ee.tex")
latexfile_hairy_cohom_eo = os.path.join(latexdir, "hairy_cohom_eo.tex")
latexfile_hairy_cohom_oe = os.path.join(latexdir, "hairy_cohom_oe.tex")
latexfile_hairy_cohom_oo = os.path.join(latexdir, "hairy_cohom_oo.tex")

latexfile_forested_vs = os.path.join(latexdir, "forested_vs.tex")
latexfile_forested_ops = os.path.join(latexdir, "forested_ops.tex")
latexfile_forested_cohom = os.path.join(latexdir, "forested_cohom.tex")

latexfile_chairy_vs = os.path.join(latexdir, "chairy_vs.tex")
latexfile_chairy_ops = os.path.join(latexdir, "chairy_ops.tex")
latexfile_chairy_cohom = os.path.join(latexdir, "chairy_cohom.tex")
latexfile_chairy_cohom_e = os.path.join(latexdir, "chairy_cohom_e.tex")
latexfile_chairy_cohom_o = os.path.join(latexdir, "chairy_cohom_o.tex")

latexfile_bichairy_vs = os.path.join(latexdir, "bichairy_vs.tex")
latexfile_bichairy_ops = os.path.join(latexdir, "bichairy_ops.tex")
latexfile_bichairy_cohom = os.path.join(latexdir, "bichairy_cohom.tex")

latexfile_forested_top_vs = os.path.join(latexdir, "forested_top_vs.tex")
latexfile_forested_top_ops = os.path.join(latexdir, "forested_top_ops.tex")
latexfile_forested_top_cohom = os.path.join(latexdir, "forested_top_cohom.tex")

latexfile_forested_nobl_top_vs = os.path.join(
    latexdir, "forested_top_vs_nobl.tex")
latexfile_forested_nobl_top_ops = os.path.join(
    latexdir, "forested_top_ops_nobl.tex")
latexfile_forested_nobl_top_cohom = os.path.join(
    latexdir, "forested_top_cohom_nobl.tex")


latexfile_forested_pre_vs = os.path.join(latexdir, "forested_pre_vs.tex")


latexfile_alldata = os.path.join(latexdir, "alldata_paper.tex")


StoreLoad.makedirs(latexdir)

alldata_tex = r"""
\documentclass{amsart}
%\usepackage{fullpage}
\usepackage[a4paper, margin=0.5in]{geometry}
\usepackage{hyperref}
\usepackage{graphicx}
\usepackage{diagbox}


\hypersetup{
    colorlinks=true, %set true if you want colored links
    linktoc=all,     %set to all if you want both sections and subsections linked
    linkcolor=blue,  %choose some color if you want links to stand out
}

\usepackage{color, colortbl}
\usepackage{array}
\usepackage{varwidth} %for the varwidth minipage environment

\definecolor{Gray}{gray}{0.9}
\definecolor{LightGray}{gray}{0.95}

\newcolumntype{g}{>{\columncolor{Gray}}c}
%\newcolumntype{M}{>{\begin{varwidth}{4cm}}c<{\end{varwidth}}} %M is for Maximal column
\newcolumntype{M}{V{3cm}}
\newcolumntype{D}{V{6cm}}
\newcolumntype{F}{p{.75cm}} % fixed width


\usepackage[table]{xcolor}

\begin{document}

\section{Ordinary}
\input{ordinary_cohom.tex}

\newpage

\section{Ordinary Merkulov}
\input{ordinaryme_cohom.tex}

\newpage

\section{Hairy}
\subsection{ee}

\begin{center}
\input{hairy_cohom_ee.tex}
\end{center}
\newpage

\subsection{eo}

\begin{center}
\input{hairy_cohom_eo.tex}
\end{center}

\newpage
\subsection{oe}

\begin{center}
\input{hairy_cohom_oe.tex}
\end{center}

\newpage
\subsection{oo}

\begin{center}
\input{hairy_cohom_oo.tex}
\end{center}

\newpage

\section{CHairy}
\input{chairy_cohom.tex}

\newpage

\section{Forested Top}
\input{forested_top_cohom.tex}

\end{document}
"""

# Refrence Euler characteristics
# from zivkovic-willwacher
ordinary_ec_evenedges = [0, 0, 0, 1, 1, 2, 1, 2, 2, 2, 1, 3, 1, 3, 4, 2, 2]
ordinary_ec_oddedges = [0, 0, 0, 1, 0, 1, -1, 1, 0,
                        0, -2, 1, 0, 0, -2, 0, -4, -3, -1, 8, 12, 27]
ordinary_ec = {True: ordinary_ec_evenedges, False: ordinary_ec_oddedges}

# param order: hairs, loops
r = "?"
# from payne-willwacher
wrhairy_ec_w2 = [[0, 0, 0, 0, 0, 0, 0, 0, -1, 4, -4],
                 [0, 1, 0, 1, 0, 1, -3, 1, -4, 9, -9],
                 [0, 0, 0, 2, 0, 4, -5, 2, -17, 9, -16],
                 [0, 0, -3, -2, 2, 12, 6, 4, -38, r, r, r],
                 [1, -3, 1, -1, r, r, r, r, r, r, r, r, r],
                 [-5, 16, r, r, r, r, r, r, r, r, r, r, r, r],
                 [r, r, r, r, r, r, r, r, r, r, r, r],
                 [r, r, r, r, r, r, r, r, r, r, r, r],
                 [r, r, r, r, r, r, r, r, r, r, r, r]]

wrhairy_ec_w1 = [[r for j in range(12)] for i in range(12)]

wrhairy_ec = {1: wrhairy_ec_w1, 2: wrhairy_ec_w2}

# from arone-turchin, hairs, loops
hairy_ec_odde_evenh = []


hairy_ec_evene_oddh = [[0, -1, -1, -1, -2, -1, -2, -2, -2, -1, -2, 0],
                       [1, 1, 1, 1, 2, 2, 3, 4, 5, 5, 7, 5],
                       [0, 0, -1, -1, -3, -4, -6, -10, -14, -17, -22, -25, -22],
                       [0, 1, 1, 2, 4, 7, 12, 20, 30, 45, 60, 79, 81, 83]]

# Color of cells depending on whether ccohomology must be zero
cell_color = {True: r" \cellcolor{LightGray}", False: ""}


def latex_table(header, data, scale=1, coltype="M", hlines=False):
    """Generates the latex Code for one table.

    :param header: list of header cells
    :type header: listr(string)
    :param data: list of rows
    :type data: list(list(str))
    :return: The Latex code of the table
    :rtype: string
    """
    colcount = len(header)
    s = "\\begin{tabular}{|g"
    for i in range(colcount-1):
        s = s + f"|{coltype}"
    s = s + "|}\n \\rowcolor{Gray}\n"

    # header
    s = s + "\\hline\n" + " & ".join(header) + "\\\\ \n" + "\\hline\n"

    for (i, row) in enumerate(data):
        if i > 0 and hlines:
            s = s+"\n \\hline \n"
        s = s + " & ".join(row) + "\\\\ \n"
        # s = s + "$" + "$ & $".join(row) + "$\\\\ \n"

    s = s + "\\hline\n\\end{tabular}\n\n"

    if scale != 1:
        s = f"\\scalebox{{ {scale} }}{{\n {s} \n }}"

    return s


def vs_dim_formatted(vs):
    if not vs.is_valid():
        return "-"
    if not vs.exists_basis_file():
        return "?"
    return str(vs.get_dimension())


def vs_dim_polynomial(vslist):
    """ Takes list of pairs (exponent, vs) """
    s = "$"
    for (exp, vs) in vslist:
        if not vs.is_valid():
            continue
        if len(s) > 1:
            s = s+"+"
        if not vs.exists_basis_file():
            s = s + f"\\text{{?}} t^{{ {exp} }} "
        else:
            s = s + f"{vs.get_dimension()}t^{{ {exp} }}"
    return s + " $"


def ops_formatted(op):
    if not op.is_valid():
        return "-"
    if not op.exists_matrix_file():
        return "?"
    if not op.exists_rank_file():
        return "R=?"
    # check if rank is mod p or rational
    r_str = "p"
    if op.exists_exact_rank():
        r_str = ""
    return f"{str(op.get_matrix_rank())} {r_str}"


def cohom_formatted(cohom_dict, tuple):
    if not tuple in cohom_dict:
        return "?"

    dim = cohom_dict[tuple]
    if not dim:
        return "?"
    else:
        return str(dim)


def cohom_formatted2(D1, D2, dim_bias=0, compute_iso=False):
    vs = D1.get_domain()
    if not vs.is_valid():
        return "-"
    if not vs.exists_basis_file():
        return "?"
    d = vs.get_dimension()

    r1 = 0
    r2 = 0
    if D1.is_valid():
        if D1.exists_rank_file():
            r1 = D1.get_matrix_rank()
        else:
            return "?"
    if D2.is_valid():
        if D2.exists_rank_file():
            r2 = D2.get_matrix_rank()
        else:
            return "?"

    # exact or not?
    r_str = "" if D1.exists_exact_rank() and D2.exists_exact_rank() else " p"
    cohomdim = d-r1-r2 + dim_bias

    # isotypical components
    isostr = ""
    if compute_iso and cohomdim > 0 and vs.get_n() >= 2:
        isostr = " (" + get_iso_string(D1, D2) + ")"

    return str(cohomdim) + r_str + isostr

# def get_forested_isostring(l, m, h, even_edges):
#     vs = ForestedGraphComplex.ForestedDegSlice(l,m,h,even_edges)
#     for i, p in enumerate(Partitions(h)):
#         P = vs.get_isotypical_projector()


iso_strings = {}


def get_iso_string(D1: SymmetricGraphComplex.SymmetricBiOperatorMatrix, D2: SymmetricGraphComplex.SymmetricBiOperatorMatrix):
    vs = D1.domain
    ret = []
    for i, p in enumerate(Partitions(vs.get_n())):
        D1iso = D1.restrict_to_isotypical_component(i)
        D2iso = D2.restrict_to_isotypical_component(i)
        isovs = D1iso.domain
        part_str = "s_{" + str(isovs.opP.rep_partition) + "}"
        if not isovs.opP.exists_rank_file():
            ret.append("?" + part_str)
            continue
        bias = - isovs.get_dimension() + isovs.get_iso_dimension()
        ret.append(cohom_formatted2(D1iso, D2iso, dim_bias=bias)
                   + part_str)
    iso_string = "$ " + ",".join(ret) + "$"
    # todo: make sure we have no hash collisions among the different complexes
    iso_strings[D1.domain] = iso_string
    return iso_string


def cohom_formatted_forested_top(D1, D2, Dc2, use_Dc2_rank=None, iso_dict=None):
    vs = D1.get_domain()
    if not vs.is_valid():
        return "-"
    if not vs.exists_basis_file():
        return "?"
    d = vs.get_dimension()

    r1 = 0
    r2 = 0
    rc2 = 0

    if Dc2.is_valid():
        if use_Dc2_rank is not None:
            if use_Dc2_rank == "?":
                return "?"
            else:
                rc2 = int(use_Dc2_rank)
        elif Dc2.exists_rank_file():
            rc2 = Dc2.get_matrix_rank()
        else:
            return "?"

    if D1.is_valid():
        if D1.exists_rank_file():
            r1 = D1.get_matrix_rank()
        else:
            return "?"
    if D2.is_valid():
        if D2.exists_rank_file():
            r2 = D2.get_matrix_rank()
        else:
            return "?"

    # exact or not?
    r_str = "" if D1.exists_exact_rank() and D2.exists_exact_rank(
    ) and Dc2.exists_exact_rank() else " p"

    # iso string
    cohomdim = d+rc2-r1-r2
    iso_str = ""
    if cohomdim > 0 and iso_dict is not None:
        fullvs = ForestedGraphComplex.ForestedDegSlice(
            vs.n_loops, vs.n_marked_edges, vs.n_hairs, vs.even_edges)
        if fullvs in iso_dict:
            iso_str = " ("+iso_dict[fullvs]+")"

    return str(cohomdim) + r_str + iso_str


def cohom_formatted_merkulov(D1, D2, Dc2):
    """The formula is the same as for the forested top complex, so we just reuse the other function."""
    return cohom_formatted_forested_top(D1, D2, Dc2)


def forested_contract_euler_rank(l, m, h, even_edges):
    """Since the contract cohomology is concentrated in top degree
    (...all vertices trivalent) the cohomology of the top
    piece of the contract operator can be found using dim counting."""
    top_v = 2*l-2+h
    ec = 0
    # contract k of the marked edges
    for k in range(1, m+1):
        vs = ForestedGraphComplex.ForestedGVS(top_v - k, l, m-k, h, even_edges)
        if vs.is_valid():
            if not vs.exists_basis_file():
                return "?"
            dim = vs.get_dimension()
            ec = ec - dim * ((-1)**k)
    return str(ec)


def eulerize(data, sign_shift=0):
    """Takes a vector of formatted dimensions (as produced by vs_dim_formatted)
    and appends an euler characteristic. """
    euler = 0
    for (i, s1) in enumerate(data):
        s = s1.split(" ")[0]
        if s == "?":
            return data + ["?"]
        elif s != "-":
            euler = euler + ((-1) ** (i+sign_shift)) * int(s)

    return data + [str(euler)]


def create_wrhairy_vs_table(v_range, l_range, h_range, w_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range] + [r"$\chi$", r"$\chi_{ref}$"]

    for w in w_range:
        for h in h_range:
            s = s + f"\n\n\\smallskip\n{w} omegas, {h} hairs \n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + eulerize(
                        [vs_dim_formatted(WRHairyGraphComplex.WRHairyGraphVS(v, l, h, w)) for v in v_range])
                    + [str(wrhairy_ec[w][h][l])])
            s = s+latex_table(header, data)
    return s


def create_wrhairy_ops_table(v_range, l_range, h_range, w_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]

    for w in w_range:
        for h in h_range:
            s = s + f"\n\n\\smallskip\n{w} omegas, {h} hairs \n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [ops_formatted(WRHairyGraphComplex.ContractEdgesGO.generate_operator(v, l, h, w)) for v in v_range])
            s = s+latex_table(header, data)
    return s

def is_wrhairy_zero(v,l,h,w):
    if w!=2:
         return False
    deg = l+v +1
    dim = 6*l - 6 + 2*h
    if l == 0:
        vcd = h-3 
        return h < 2 or deg < dim - vcd or deg > dim
    elif h == 0:
        vcd = 4*l-5
        return deg < dim - vcd + 1 or deg > dim
    else:
        vcd = 4*l-4+h
        return deg < dim - vcd +(1 if l==1 else 0)  or deg > dim 

def create_wrhairy_cohom_table(v_range, l_range, h_range, w_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]

    for w in w_range:
        for h in h_range:
            s = s + f"\n\\smallskip\n{w} omegas, {h} hairs \n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [cohom_formatted2(
                        WRHairyGraphComplex.ContractEdgesGO.generate_operator(
                            v, l, h, w),
                        WRHairyGraphComplex.ContractEdgesGO.generate_operator(v+1, l, h, w))
                        + cell_color[is_wrhairy_zero(v, l, h, w)] for v in v_range])
            s = s+latex_table(header, data)
    return s


def is_ordinary_zero(v, l):
    return (v < l+1) or (v > 2*l-2)


def create_ordinary_vs_table(v_range, l_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range] + [r"$\chi$", r"$\chi_{ref}$"]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges \n\n"
        data = []
        for l in l_range:
            ref_ec = ordinary_ec[even_edges][l]
            if not even_edges:
                ref_ec = ((-1)**(l+1)) * ref_ec
            data.append(
                [str(l)] + eulerize(
                    [vs_dim_formatted(OrdinaryGraphComplex.OrdinaryGVS(
                        v, l, even_edges)) + cell_color[is_ordinary_zero(v, l)]
                        for v in v_range]
                ) + [str(ref_ec)]
            )
        s = s+latex_table(header, data, scale=0.75)
    return s


def create_ordinary_ops_table(v_range, l_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges \n\n"
        data = []
        for l in l_range:
            data.append(
                [str(l)] + [ops_formatted(OrdinaryGraphComplex.ContractEdgesGO.generate_operator(v, l, even_edges)) for v in v_range])
        s = s+latex_table(header, data)
    return s


def create_ordinary_cohom_table(v_range, l_range):
    s = ""

    header = ["\\diagbox{l}{v}"] + \
        [str(v) for v in v_range] + [r"$\chi$", r"$\chi_{ref}$"]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges \n\n"
        data = []
        for l in l_range:
            ref_ec = ordinary_ec[even_edges][l]
            if not even_edges:
                ref_ec = ((-1)**(l+1)) * ref_ec
            data.append(
                [str(l)] + eulerize(
                    [cohom_formatted2(
                        OrdinaryGraphComplex.ContractEdgesGO.generate_operator(
                            v, l, even_edges),
                        OrdinaryGraphComplex.ContractEdgesGO.generate_operator(
                            v+1, l, even_edges)
                    ) + cell_color[is_ordinary_zero(v, l)] for v in v_range])
                + [str(ref_ec)]
            )
        s = s+latex_table(header, data)
    return s


def create_ordinaryme_vs_table(v_range, l_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges \n\n"
        data = []
        for l in l_range:
            data.append(
                [str(l)] +
                [vs_dim_formatted(OrdinaryMerkulovComplex.OrdinaryMerkulovGVS(
                    v, l, even_edges, 3456)) + cell_color[is_ordinary_zero(v, l)]
                 for v in v_range]
            )
        s = s+latex_table(header, data, scale=0.75)
    return s


def create_ordinaryme_ops_table(v_range, l_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges \n\n"
        data = []
        for l in l_range:
            data.append(
                [str(l)] + [ops_formatted(OrdinaryMerkulovComplex.ContractEdgesGO.generate_operator(v, l, even_edges)) for v in v_range])
        s = s+latex_table(header, data)
    return s


def create_ordinaryme_cohom_table(v_range, l_range):
    s = ""

    header = ["\\diagbox{l}{v}"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges \n\n"
        data = []
        for l in l_range:
            data.append(
                [str(l)] + [cohom_formatted_merkulov(
                    OrdinaryMerkulovComplex.ContractEdgesGO.generate_operator(
                        v, l, even_edges),
                    OrdinaryMerkulovComplex.ContractEdgesGO.generate_operator(
                        v+1, l, even_edges),
                    OrdinaryMerkulovComplex.ContractEdgesGO.generate_operator(
                        v+1, l, even_edges, False)
                ) + cell_color[is_ordinary_zero(v, l)] for v in v_range])
        s = s+latex_table(header, data)
    return s


def is_hairy_zero(v, l, h):
    """Determines whether the cohomology in the given degree is zero by abstract reasons."""
    return (v < l+h-2) or (v > 2*l-2+h)


def create_hairy_cohom_table(v_range, hl_pairs, even_edges, even_hairs):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    # for even_edges in [True, False]:
    #     for even_hairs in [True, False]:
    # s = s + "\n\n\\smallskip\n\n" #+ \
        # ("even" if even_edges else "odd") + " edges, " + \
        # ("even" if even_hairs else "odd") + " hairs \n\n"
    for h, l_range in hl_pairs:
        s = s + f"\n\n\\smallskip\n\n{h} hairs\n\n"
        data = []
        for l in l_range:
            data.append(
                [str(l)] + [cohom_formatted2(
                    HairyGraphComplex.ContractEdgesGO.generate_operator(
                        v, l, h, even_edges, even_hairs),
                    HairyGraphComplex.ContractEdgesGO.generate_operator(
                        v+1, l, h, even_edges, even_hairs)
                ) + cell_color[is_hairy_zero(v, l, h)] for v in v_range])
        s = s+latex_table(header, data,scale=0.75, coltype="F")
    return s


def create_hairy_vs_table(v_range, l_range, h_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        for even_hairs in [True, False]:
            s = s + "\n\n\\smallskip\n" + \
                ("even" if even_edges else "odd") + " edges, " + \
                ("even" if even_hairs else "odd") + " hairs \n\n"
            for h in h_range:
                s = s + f"\n{h} hairs\n\n"
                data = []
                for l in l_range:
                    data.append(
                        [str(l)] + [vs_dim_formatted(
                            HairyGraphComplex.HairyGraphVS(
                                v, l, h, even_edges, even_hairs)
                        ) + cell_color[is_hairy_zero(v, l, h)] for v in v_range])
                s = s+latex_table(header, data)
    return s


def create_hairy_ops_table(v_range, l_range, h_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        for even_hairs in [True, False]:
            s = s + "\n\n\\smallskip\n" + \
                ("even" if even_edges else "odd") + " edges, " + \
                ("even" if even_hairs else "odd") + " hairs \n\n"
            for h in h_range:
                s = s + f"\n{h} hairs\n\n"
                data = []
                for l in l_range:
                    data.append(
                        [str(l)] + [ops_formatted(
                            HairyGraphComplex.ContractEdgesGO.generate_operator(
                                v, l, h, even_edges, even_hairs)
                        ) for v in v_range])
                s = s+latex_table(header, data)
    return s


def create_bichairy_vs_table(v_range, l_range, h_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        for even_hairs_a in [True, False]:
            for even_hairs_b in [True, False]:
                s = s + "\n\n\\smallskip\n" + \
                    ("even" if even_edges else "odd") + " edges, " + \
                    ("even" if even_hairs_a else "odd") + " hairs a " + \
                    ("even" if even_hairs_b else "odd") + " hairs b \n\n"
                for h in h_range:
                    for ha in range(h):
                        s = s + f"\n{h} hairs ({ha}+{h-ha})\n\n"
                        data = []
                        for l in l_range:
                            data.append(
                                [str(l)] + [vs_dim_formatted(
                                    BiColoredHairyGraphComplex.BiColoredHairyGraphVS(
                                        v, l, ha, h-ha, even_edges, even_hairs_a, even_hairs_b)
                                ) for v in v_range])
                        s = s+latex_table(header, data)
    return s


def create_bichairy_ops_table(v_range, l_range, h_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        for even_hairs_a in [True, False]:
            for even_hairs_b in [True, False]:
                s = s + "\n\n\\smallskip\n" + \
                    ("even" if even_edges else "odd") + " edges, " + \
                    ("even" if even_hairs_a else "odd") + " hairs a " + \
                    ("even" if even_hairs_b else "odd") + " hairs b \n\n"
                for h in h_range:
                    for ha in range(h):
                        s = s + f"\n\n{h} hairs ({ha}+{h-ha})\n\n"
                        data = []
                        for l in l_range:
                            data.append(
                                [str(l)] + [ops_formatted(
                                    BiColoredHairyGraphComplex.ContractEdgesGO.generate_operator(
                                        v, l, ha, h-ha, even_edges, even_hairs_a, even_hairs_b)
                                ) for v in v_range])
                        s = s+latex_table(header, data)
    return s


def create_bichairy_cohom_table(v_range, l_range, h_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        for even_hairs_a in [True, False]:
            for even_hairs_b in [True, False]:
                s = s + "\n\n\\smallskip\n" + \
                    ("even" if even_edges else "odd") + " edges, " + \
                    ("even" if even_hairs_a else "odd") + " hairs a " + \
                    ("even" if even_hairs_b else "odd") + " hairs b \n\n"
                for h in h_range:
                    for ha in range(h):
                        s = s + f"\n\n{h} hairs ({ha}+{h-ha})\n\n"
                        data = []
                        for l in l_range:
                            data.append(
                                [str(l)] + [cohom_formatted2(
                                    BiColoredHairyGraphComplex.ContractEdgesGO.generate_operator(
                                        v, l, ha, h-ha, even_edges, even_hairs_a, even_hairs_b),
                                    BiColoredHairyGraphComplex.ContractEdgesGO.generate_operator(
                                        v+1, l, ha, h-ha, even_edges, even_hairs_a, even_hairs_b)
                                ) for v in v_range])
                        s = s+latex_table(header, data)
    return s


def create_chairy_vs_table(v_range, l_range, h_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n{h} hairs\n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [vs_dim_formatted(
                        CHairyGraphComplex.CHairyGraphVS(
                            v, l, h, even_edges)
                    ) for v in v_range])
            s = s+latex_table(header, data)
    return s


def create_chairy_ops_table(v_range, l_range, h_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n\n{h} hairs\n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [ops_formatted(
                        CHairyGraphComplex.ContractEdgesGO.generate_operator(
                            v, l, h, even_edges)
                    ) for v in v_range])
            s = s+latex_table(header, data)
    return s


def create_chairy_cohom_table(v_range, l_range, h_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n\n{h} hairs\n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [cohom_formatted2(
                        CHairyGraphComplex.ContractEdgesGO.generate_operator(
                            v, l, h, even_edges),
                        CHairyGraphComplex.ContractEdgesGO.generate_operator(
                            v+1, l, h, even_edges),
                        compute_iso=True
                    ) for v in v_range])
            s = s+latex_table(header, data)
    return s


def create_forested_pre_vs_table(v_range, l_range, m_range, h_range):
    s = ""

    header = ["v,l"] + [str(l) for l in l_range]
    for h in h_range:
        s = s + f"\n\n{h} hairs\n\n"
        data = []
        for v in v_range:
            data.append(
                [str(v)] + [vs_dim_polynomial(
                    [
                        (m, ForestedGraphComplex.PreForestedGVS(
                            v, l, m, h))
                        for m in m_range]
                ) for l in l_range])
        s = s+latex_table(header, data, scale=.5, coltype="D", hlines=True)
    return s


def create_forested_vs_table(l_range, m_range, h_range):
    s = ""

    header = ["l,m"] + [str(m) for m in m_range]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n\n{h} hairs\n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [vs_dim_formatted(
                        ForestedGraphComplex.ForestedDegSlice(
                            l, m, h, even_edges)
                    ) for m in m_range])
            s = s+latex_table(header, data)
    return s


def create_forested_ops_table(l_range, m_range, h_range):
    s = ""

    header = ["l,m"] + [str(m) for m in m_range]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n\n{h} hairs\n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [ops_formatted(
                        ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(
                            l, m, h, even_edges)
                    ) for m in m_range])
            s = s+latex_table(header, data)
    return s


def create_forested_cohom_table(l_range, m_range, h_range):
    s = ""

    header = ["l,m"] + [str(m) for m in m_range]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n\n{h} hairs\n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [cohom_formatted2(
                        ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(
                            l, m, h, even_edges),
                        ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(
                            l, m+1, h, even_edges),
                        compute_iso=True
                    ) for m in m_range])
            s = s+latex_table(header, data)
    return s


def create_forested_top_vs_table(l_range, m_range, h_range):
    s = ""

    header = ["l,m"] + [str(m) for m in m_range]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            for topn in [1, 2]:
                s = s + f"\n\n{h} hairs, {topn} topn\n\n"
                data = []
                for l in l_range:
                    data.append(
                        [str(l)] + [vs_dim_formatted(
                            ForestedGraphComplex.ForestedTopDegSlice(
                                l, m, h, even_edges, topn)
                        ) for m in m_range])
                s = s+latex_table(header, data)
    return s


def create_forested_top_ops_table(l_range, m_range, h_range):
    s = ""

    header = ["l,m"] + [str(m) for m in m_range]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n\n{h} hairs, contractunmarktop \n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [ops_formatted(
                        ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(
                            l, m, h, even_edges)
                    ) for m in m_range])
            s = s+latex_table(header, data)

            s = s + f"\n{h} hairs, contract \n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [ops_formatted(
                        ForestedGraphComplex.ContractEdgesGO.generate_operator(
                            2*l-2+h, l, m, h, even_edges)
                    ) + " (" +
                        forested_contract_euler_rank(l, m, h, even_edges)
                        + ")" for m in m_range])
            s = s+latex_table(header, data, scale=0.75)

    return s


def create_forested_top_cohom_table(l_range, m_range, h_range):
    s = ""

    header = ["l,m"] + [str(m) for m in m_range]
    for even_edges in [True, False]:
        s = s + "\n\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n\n{h} hairs, contractunmarktop \n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [cohom_formatted_forested_top(
                        ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(
                            l, m, h, even_edges),
                        ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(
                            l, m+1, h, even_edges),
                        ForestedGraphComplex.ContractEdgesGO.generate_operator(
                            2*l-2+h, l, m+1, h, even_edges),
                        use_Dc2_rank=forested_contract_euler_rank(
                            l, m+1, h, even_edges),
                        iso_dict=iso_strings
                    ) for m in m_range])
            s = s+latex_table(header, data)

    return s


def write_tables():
    # Generate tables
    print("Ordinary....")
    # s = create_ordinary_vs_table(range(4, 24), range(3, 13))
    # with open(latexfile_ordinary_vs, 'w') as f:
    #     f.write(s)

    # s = create_ordinary_ops_table(range(4, 24), range(3, 13))
    # with open(latexfile_ordinary_ops, 'w') as f:
    #     f.write(s)

    # s = create_ordinary_cohom_table(range(4, 22), range(3, 12))
    # with open(latexfile_ordinary_cohom, 'w') as f:
    #     f.write(s)

    print("Ordinary Merkulov....")
    # s = create_ordinaryme_vs_table(range(4, 24), range(3, 13))
    # with open(latexfile_ordinaryme_vs, 'w') as f:
    #     f.write(s)

    # s = create_ordinaryme_ops_table(range(4, 25), range(3, 13))
    # with open(latexfile_ordinaryme_ops, 'w') as f:
    #     f.write(s)

    # s = create_ordinaryme_cohom_table(range(4, 22), range(3, 12))
    # with open(latexfile_ordinaryme_cohom, 'w') as f:
    #     f.write(s)

    print("Hairy....")
    # s = create_hairy_vs_table(range(22), range(12), range(1, 9))
    # with open(latexfile_hairy_vs, 'w') as f:
    #     f.write(s)

    # s = create_hairy_ops_table(range(20), range(12), range(1, 9))
    # with open(latexfile_hairy_ops, 'w') as f:
    #     f.write(s)
    hl_pairs = [(1,range(10)), (2,range(9)),(3,range(8)),(4,range(7)),(5,range(7)),(6,range(6)),(7,range(6)),(8,range(6))]
    s = create_hairy_cohom_table(range(18), hl_pairs, True, True)
    with open(latexfile_hairy_cohom_ee, 'w') as f:
        f.write(s)
    s = create_hairy_cohom_table(range(18), hl_pairs, True, False)
    with open(latexfile_hairy_cohom_eo, 'w') as f:
        f.write(s)
    s = create_hairy_cohom_table(range(18), hl_pairs, False, True)
    with open(latexfile_hairy_cohom_oe, 'w') as f:
        f.write(s)
    s = create_hairy_cohom_table(range(18), hl_pairs, False, False)
    with open(latexfile_hairy_cohom_oo, 'w') as f:
        f.write(s)

    print("CHairy....")
    # s = create_chairy_vs_table(range(20), range(12), range(6))
    # with open(latexfile_chairy_vs, 'w') as f:
    #     f.write(s)

    # s = create_chairy_ops_table(range(20), range(12), range(6))
    # with open(latexfile_chairy_ops, 'w') as f:
    #     f.write(s)

    s = create_chairy_cohom_table(range(20), range(12), range(6))
    with open(latexfile_chairy_cohom, 'w') as f:
        f.write(s)

    # print("BiColoredHairy....")
    # s = create_bichairy_vs_table(range(25), range(12), range(6))
    # with open(latexfile_bichairy_vs, 'w') as f:
    #     f.write(s)

    # s = create_bichairy_ops_table(range(25), range(12), range(6))
    # with open(latexfile_bichairy_ops, 'w') as f:
    #     f.write(s)

    # s = create_bichairy_cohom_table(range(25), range(12), range(6))
    # with open(latexfile_bichairy_cohom, 'w') as f:
    #     f.write(s)

    # print("WRHairy....")
    # s = create_wrhairy_vs_table(range(25), range(11), range(8), range(1, 3))
    # with open(latexfile_wrhairy_vs, 'w') as f:
    #     f.write(s)

    # s = create_wrhairy_ops_table(range(25), range(11), range(8), range(1, 3))
    # with open(latexfile_wrhairy_ops, 'w') as f:
    #     f.write(s)

    # s = create_wrhairy_cohom_table(range(21), range(11), range(8), range(1, 3))
    # with open(latexfile_wrhairy_cohom, 'w') as f:
    #     f.write(s)

    # print("Forested....")
    # s = create_forested_vs_table(range(9), range(20), range(6))
    # with open(latexfile_forested_vs, 'w') as f:
    #     f.write(s)

    # s = create_forested_pre_vs_table(range(20), range(9), range(20), range(6))
    # with open(latexfile_forested_pre_vs, 'w') as f:
    #     f.write(s)

    # s = create_forested_ops_table(range(9), range(20), range(6))
    # with open(latexfile_forested_ops, 'w') as f:
    #     f.write(s)

    # s = create_forested_cohom_table(range(9), range(20), range(6))
    # with open(latexfile_forested_cohom, 'w') as f:
    #     f.write(s)

    print("Forested Top BL....")
    ForestedGraphComplex.use_bridgeless = True
    ForestedGraphComplex.graph_type = "forestedbl"
    # s = create_forested_top_vs_table(range(9), range(20), range(6))
    # with open(latexfile_forested_top_vs, 'w') as f:
    #     f.write(s)

    # s = create_forested_top_ops_table(range(9), range(20), range(6))
    # with open(latexfile_forested_top_ops, 'w') as f:
    #     f.write(s)

    # s = create_forested_top_cohom_table(range(9), range(20), range(6))
    # with open(latexfile_forested_top_cohom, 'w') as f:
    #     f.write(s)

    # print("Forested Top noBL....")
    # ForestedGraphComplex.use_bridgeless = False
    # ForestedGraphComplex.graph_type = "forested"
    # s = create_forested_top_vs_table(range(9), range(20), range(6))
    # with open(latexfile_forested_nobl_top_vs, 'w') as f:
    #     f.write(s)

    # s = create_forested_top_ops_table(range(9), range(20), range(6))
    # with open(latexfile_forested_nobl_top_ops, 'w') as f:
    #     f.write(s)

    # s = create_forested_top_cohom_table(range(9), range(20), range(6))
    # with open(latexfile_forested_nobl_top_cohom, 'w') as f:
    #     f.write(s)


def write_alldata():
    with open(latexfile_alldata, 'w') as f:
        f.write(alldata_tex)


def all_export():
    write_tables()
    write_alldata()


all_export()


#  def print_dim_and_eulerchar(self):
#         for w in self.w_range:
#             for h in self.h_range:
#                 for l in self.l_range:
#                     ds = [WRHairyGraphVS(v, l, h, w).get_dimension()
#                           for v in self.v_range]
#                     eul = sum([(1 if j % 2 == 0 else -1) *
#                               d for j, d in enumerate(ds)])
#                     print("Dimensions (w,h,l) ", w,
#                           h, l, ":", ds, "Euler", eul)

#     def print_cohomology_dim(self):
#         for w in self.w_range:
#             for h in self.h_range:
#                 for l in self.l_range:
#                     cohomdict = {}
#                     for v in self.v_range:
#                         D1 = ContractEdgesGO.generate_operator(v, l, h, w)
#                         D2 = ContractEdgesGO.generate_operator(v+1, l, h, w)
#                         try:
#                             d = WRHairyGraphVS(v, l, h, w).get_dimension()
#                             r1 = D1.get_matrix_rank()
#                             r2 = D2.get_matrix_rank()
#                             cohomdict[v] = d-r1-r2
#                         except:
#                             pass

#                     print("Cohomology Dimensions (w,h,l) ",
#                           w, h, l, ":", cohomdict)
