"""This script generates the latex tables for the paper."""

import OrdinaryGraphComplex
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

latexdir = os.path.join(Parameters.plots_dir, "latex")

latexfile_wrhairy_vs = os.path.join(latexdir, "wrhairy_vs.tex")
latexfile_wrhairy_ops = os.path.join(latexdir, "wrhairy_ops.tex")
latexfile_wrhairy_cohom = os.path.join(latexdir, "wrhairy_cohom.tex")

latexfile_ordinary_vs = os.path.join(latexdir, "ordinary_vs.tex")
latexfile_ordinary_ops = os.path.join(latexdir, "ordinary_ops.tex")
latexfile_ordinary_cohom = os.path.join(latexdir, "ordinary_cohom.tex")

latexfile_hairy_vs = os.path.join(latexdir, "hairy_vs.tex")
latexfile_hairy_ops = os.path.join(latexdir, "hairy_ops.tex")
latexfile_hairy_cohom = os.path.join(latexdir, "hairy_cohom.tex")

latexfile_forested_vs = os.path.join(latexdir, "forested_vs.tex")
latexfile_forested_ops = os.path.join(latexdir, "forested_ops.tex")
latexfile_forested_cohom = os.path.join(latexdir, "forested_cohom.tex")

latexfile_chairy_vs = os.path.join(latexdir, "chairy_vs.tex")
latexfile_chairy_ops = os.path.join(latexdir, "chairy_ops.tex")
latexfile_chairy_cohom = os.path.join(latexdir, "chairy_cohom.tex")

latexfile_bichairy_vs = os.path.join(latexdir, "bichairy_vs.tex")
latexfile_bichairy_ops = os.path.join(latexdir, "bichairy_ops.tex")
latexfile_bichairy_cohom = os.path.join(latexdir, "bichairy_cohom.tex")

latexfile_forested_top_vs = os.path.join(latexdir, "forested_top_vs.tex")
latexfile_forested_top_ops = os.path.join(latexdir, "forested_top_ops.tex")
latexfile_forested_top_cohom = os.path.join(latexdir, "forested_top_cohom.tex")

latexfile_forested_pre_vs = os.path.join(latexdir, "forested_pre_vs.tex")


latexfile_alldata = os.path.join(latexdir, "alldata.tex")


StoreLoad.makedirs(latexdir)

alldata_tex = r"""
\documentclass{amsart}
%\usepackage{fullpage}
\usepackage[a4paper, landscape, margin=0.5in]{geometry}
\usepackage{hyperref}
\usepackage{graphicx}


\hypersetup{
    colorlinks=true, %set true if you want colored links
    linktoc=all,     %set to all if you want both sections and subsections linked
    linkcolor=blue,  %choose some color if you want links to stand out
}

\usepackage{color, colortbl}
\usepackage{array}
\usepackage{varwidth} %for the varwidth minipage environment

\definecolor{Gray}{gray}{0.9}

\newcolumntype{g}{>{\columncolor{Gray}}c}
%\newcolumntype{M}{>{\begin{varwidth}{4cm}}c<{\end{varwidth}}} %M is for Maximal column
\newcolumntype{M}{V{3cm}}
\newcolumntype{D}{V{6cm}}

\begin{document}

\section{Ordinary}

\subsection{VS Dimensions}
\input{ordinary_vs.tex}
 
\subsection{Operator ranks}
\input{ordinary_ops.tex}

\subsection{Cohomology}
\input{ordinary_cohom.tex}

\newpage

\section{Hairy}

\subsection{VS Dimensions}
\input{hairy_vs.tex}
 
\subsection{Operator ranks}
\input{hairy_ops.tex}

\subsection{Cohomology}
\input{hairy_cohom.tex}

\newpage

\section{CHairy}

\subsection{VS Dimensions}
\input{chairy_vs.tex}
 
\subsection{Operator ranks}
\input{chairy_ops.tex}

\subsection{Cohomology}
\input{chairy_cohom.tex}

\newpage

\section{BiCHairy}

\subsection{VS Dimensions}
\input{bichairy_vs.tex}
 
\subsection{Operator ranks}
\input{bichairy_ops.tex}

\subsection{Cohomology}
\input{bichairy_cohom.tex}

\newpage

\section{WRHairy}

\subsection{VS Dimensions}
\input{wrhairy_vs.tex}
 
\subsection{Operator ranks}
\input{wrhairy_ops.tex}

\subsection{Cohomology}
\input{wrhairy_cohom.tex}

\newpage

\section{Forested}

\subsection{PreVS Dimensions}
\input{forested_pre_vs.tex}

\subsection{VS Dimensions}
\input{forested_vs.tex}
 
\subsection{Operator ranks}
\input{forested_ops.tex}

\subsection{Cohomology}
\input{forested_cohom.tex}

\newpage

\section{Forested Top}

\subsection{VS Dimensions}
\input{forested_top_vs.tex}
 
\subsection{Operator ranks}
\input{forested_top_ops.tex}

\subsection{Cohomology}
\input{forested_top_cohom.tex}

\end{document}
"""


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


def wrhairy_vs_dim_poly(l, h, w):
    max_vertices = 2*l-2+h


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
        if len(s) >1:
            s=s+"+"
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


def cohom_formatted2(D1, D2):
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

    return str(d-r1-r2) + r_str

def cohom_formatted_forested_top(D1, D2, Dc2):
    vs = D1.get_domain()
    if not vs.is_valid():
        return "-"
    if not vs.exists_basis_file():
        return "?"
    d = vs.get_dimension()

    r1 = 0
    r2 = 0
    rc2 = 0
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

    if Dc2.is_valid():
        if Dc2.exists_rank_file():
            rc2 = Dc2.get_matrix_rank()
        else:
            return "?"

    # exact or not?
    r_str = "" if D1.exists_exact_rank() and D2.exists_exact_rank() and Dc2.exists_exact_rank() else " p"

    return str(d+rc2-r1-r2) + r_str


def create_wrhairy_vs_table(v_range, l_range, h_range, w_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]

    for w in w_range:
        for h in h_range:
            s = s + f"\n\\smallskip\n{w} omegas, {h} hairs \n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [vs_dim_formatted(WRHairyGraphComplex.WRHairyGraphVS(v, l, h, w)) for v in v_range])
            s = s+latex_table(header, data)
    return s


def create_wrhairy_ops_table(v_range, l_range, h_range, w_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]

    for w in w_range:
        for h in h_range:
            s = s + f"\n\\smallskip\n{w} omegas, {h} hairs \n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [ops_formatted(WRHairyGraphComplex.ContractEdgesGO.generate_operator(v, l, h, w)) for v in v_range])
            s = s+latex_table(header, data)
    return s


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
                        for v in v_range])
            s = s+latex_table(header, data)
    return s


def create_ordinary_vs_table(v_range, l_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        s = s + "\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges \n\n"
        data = []
        for l in l_range:
            data.append(
                [str(l)] + [vs_dim_formatted(OrdinaryGraphComplex.OrdinaryGVS(v, l, even_edges)) for v in v_range])
        s = s+latex_table(header, data)
    return s


def create_ordinary_ops_table(v_range, l_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        s = s + "\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges \n\n"
        data = []
        for l in l_range:
            data.append(
                [str(l)] + [ops_formatted(OrdinaryGraphComplex.ContractEdgesGO.generate_operator(v, l, even_edges)) for v in v_range])
        s = s+latex_table(header, data)
    return s


def create_ordinary_cohom_table(v_range, l_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        s = s + "\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges \n\n"
        data = []
        for l in l_range:
            data.append(
                [str(l)] + [cohom_formatted2(
                    OrdinaryGraphComplex.ContractEdgesGO.generate_operator(
                        v, l, even_edges),
                    OrdinaryGraphComplex.ContractEdgesGO.generate_operator(
                        v+1, l, even_edges)
                ) for v in v_range])
        s = s+latex_table(header, data)
    return s


def create_hairy_cohom_table(v_range, l_range, h_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        for even_hairs in [True, False]:
            s = s + "\n\\smallskip\n" + \
                ("even" if even_edges else "odd") + " edges, " + \
                ("even" if even_hairs else "odd") + " hairs \n\n"
            for h in h_range:
                s = s + f"\n{h} hairs\n\n"
                data = []
                for l in l_range:
                    data.append(
                        [str(l)] + [cohom_formatted2(
                            HairyGraphComplex.ContractEdgesGO.generate_operator(
                                v, l, h, even_edges, even_hairs),
                            HairyGraphComplex.ContractEdgesGO.generate_operator(
                                v+1, l, h, even_edges, even_hairs)
                        ) for v in v_range])
                s = s+latex_table(header, data)
    return s


def create_hairy_vs_table(v_range, l_range, h_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        for even_hairs in [True, False]:
            s = s + "\n\\smallskip\n" + \
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
                        ) for v in v_range])
                s = s+latex_table(header, data)
    return s


def create_hairy_ops_table(v_range, l_range, h_range):
    s = ""

    header = ["l,v"] + [str(v) for v in v_range]
    for even_edges in [True, False]:
        for even_hairs in [True, False]:
            s = s + "\n\\smallskip\n" + \
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
                s = s + "\n\\smallskip\n" + \
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
                s = s + "\n\\smallskip\n" + \
                    ("even" if even_edges else "odd") + " edges, " + \
                    ("even" if even_hairs_a else "odd") + " hairs a " + \
                    ("even" if even_hairs_b else "odd") + " hairs b \n\n"
                for h in h_range:
                    for ha in range(h):
                        s = s + f"\n{h} hairs ({ha}+{h-ha})\n\n"
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
                s = s + "\n\\smallskip\n" + \
                    ("even" if even_edges else "odd") + " edges, " + \
                    ("even" if even_hairs_a else "odd") + " hairs a " + \
                    ("even" if even_hairs_b else "odd") + " hairs b \n\n"
                for h in h_range:
                    for ha in range(h):
                        s = s + f"\n{h} hairs ({ha}+{h-ha})\n\n"
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
        s = s + "\n\\smallskip\n" + \
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
        s = s + "\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n{h} hairs\n\n"
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
        s = s + "\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n{h} hairs\n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [cohom_formatted2(
                        CHairyGraphComplex.ContractEdgesGO.generate_operator(
                            v, l, h, even_edges),
                        CHairyGraphComplex.ContractEdgesGO.generate_operator(
                            v+1, l, h, even_edges)
                    ) for v in v_range])
            s = s+latex_table(header, data)
    return s

def create_forested_pre_vs_table(v_range, l_range, m_range, h_range):
    s = ""

    header = ["v,l"] + [str(l) for l in l_range]
    for h in h_range:
        s = s + f"\n{h} hairs\n\n"
        data = []
        for v in v_range:
            data.append(
                [str(v)] + [vs_dim_polynomial(
                    [
                    (m, ForestedGraphComplex.PreForestedGVS(
                        v, l, m, h))
                    for m in m_range ]
                ) for l in l_range])
        s = s+latex_table(header, data, scale=.5, coltype="D", hlines=True)
    return s

def create_forested_vs_table(l_range, m_range, h_range):
    s = ""

    header = ["l,m"] + [str(m) for m in m_range]
    for even_edges in [True, False]:
        s = s + "\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n{h} hairs\n\n"
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
        s = s + "\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n{h} hairs\n\n"
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
        s = s + "\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n{h} hairs\n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [cohom_formatted2(
                        ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(
                            l, m, h, even_edges),
                        ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(
                            l, m+1, h, even_edges)
                    ) for m in m_range])
            s = s+latex_table(header, data)
    return s

def create_forested_top_vs_table(l_range, m_range, h_range):
    s = ""

    header = ["l,m"] + [str(m) for m in m_range]
    for even_edges in [True, False]:
        s = s + "\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            for topn in [1,2]:
                s = s + f"\n{h} hairs, {topn} topn\n\n"
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
        s = s + "\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n{h} hairs, contractunmarktop \n\n"
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
                    ) for m in m_range])
            s = s+latex_table(header, data)

    return s

def create_forested_top_cohom_table(l_range, m_range, h_range):
    s = ""

    header = ["l,m"] + [str(m) for m in m_range]
    for even_edges in [True, False]:
        s = s + "\n\\smallskip\n" + \
            ("even" if even_edges else "odd") + " edges\n\n "
        for h in h_range:
            s = s + f"\n{h} hairs, contractunmarktop \n\n"
            data = []
            for l in l_range:
                data.append(
                    [str(l)] + [cohom_formatted_forested_top(
                        ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(
                            l, m, h, even_edges),
                        ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(
                            l, m+1, h, even_edges),
                        ForestedGraphComplex.ContractEdgesGO.generate_operator(
                            2*l-2+h, l, m+1, h, even_edges)
                    ) for m in m_range])
            s = s+latex_table(header, data)

    return s

def write_tables():
    # Generate tables
    print("WRHairy....")
    s = create_wrhairy_vs_table(range(25), range(9), range(6), range(1, 3))
    with open(latexfile_wrhairy_vs, 'w') as f:
        f.write(s)

    s = create_wrhairy_ops_table(range(25), range(9), range(6), range(1, 3))
    with open(latexfile_wrhairy_ops, 'w') as f:
        f.write(s)

    s = create_wrhairy_cohom_table(range(25), range(9), range(6), range(1, 3))
    with open(latexfile_wrhairy_cohom, 'w') as f:
        f.write(s)

    print("Ordinary....")
    s = create_ordinary_vs_table(range(25), range(15))
    with open(latexfile_ordinary_vs, 'w') as f:
        f.write(s)

    s = create_ordinary_ops_table(range(25), range(15))
    with open(latexfile_ordinary_ops, 'w') as f:
        f.write(s)

    s = create_ordinary_cohom_table(range(25), range(15))
    with open(latexfile_ordinary_cohom, 'w') as f:
        f.write(s)

    print("Hairy....")
    s = create_hairy_vs_table(range(25), range(12), range(6))
    with open(latexfile_hairy_vs, 'w') as f:
        f.write(s)

    s = create_hairy_ops_table(range(25), range(12), range(6))
    with open(latexfile_hairy_ops, 'w') as f:
        f.write(s)

    s = create_hairy_cohom_table(range(25), range(12), range(6))
    with open(latexfile_hairy_cohom, 'w') as f:
        f.write(s)

    print("CHairy....")
    s = create_chairy_vs_table(range(20), range(12), range(6))
    with open(latexfile_chairy_vs, 'w') as f:
        f.write(s)

    s = create_chairy_ops_table(range(20), range(12), range(6))
    with open(latexfile_chairy_ops, 'w') as f:
        f.write(s)

    s = create_chairy_cohom_table(range(20), range(12), range(6))
    with open(latexfile_chairy_cohom, 'w') as f:
        f.write(s)

    print("BiColoredHairy....")
    s = create_bichairy_vs_table(range(25), range(12), range(6))
    with open(latexfile_bichairy_vs, 'w') as f:
        f.write(s)

    s = create_bichairy_ops_table(range(25), range(12), range(6))
    with open(latexfile_bichairy_ops, 'w') as f:
        f.write(s)

    s = create_bichairy_cohom_table(range(25), range(12), range(6))
    with open(latexfile_bichairy_cohom, 'w') as f:
        f.write(s)

    print("Forested....")
    s = create_forested_vs_table(range(9), range(20), range(6))
    with open(latexfile_forested_vs, 'w') as f:
        f.write(s)
        
    s = create_forested_pre_vs_table(range(20), range(9), range(20), range(6))
    with open(latexfile_forested_pre_vs, 'w') as f:
        f.write(s)

    s = create_forested_ops_table(range(9), range(20), range(6))
    with open(latexfile_forested_ops, 'w') as f:
        f.write(s)

    s = create_forested_cohom_table(range(9), range(20), range(6))
    with open(latexfile_forested_cohom, 'w') as f:
        f.write(s)

    print("Forested Top....")
    s = create_forested_top_vs_table(range(9), range(20), range(6))
    with open(latexfile_forested_top_vs, 'w') as f:
        f.write(s)

    s = create_forested_top_ops_table(range(9), range(20), range(6))
    with open(latexfile_forested_top_ops, 'w') as f:
        f.write(s)

    s = create_forested_top_cohom_table(range(9), range(20), range(6))
    with open(latexfile_forested_top_cohom, 'w') as f:
        f.write(s)


def write_alldata():
    with open(latexfile_alldata, 'w') as f:
        f.write(alldata_tex)


write_tables()
write_alldata()

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
