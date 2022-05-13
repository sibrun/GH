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

StoreLoad.makedirs(latexdir)

alldata_tex =r"""
\documentclass{amsart}
\usepackage{fullpage}
\usepackage{hyperref}
\begin{document}

\section{WRHairy}

\subsection{VS Dimensions}
\input{wrhairy_vs.tex}
 
\subsection{Operators}
\input{wrhairy_ops.tex}

\subsection{Cohomology}
\input{wrhairy_cohom.tex}

\end{document}
"""


def latex_table(header, data):
    """Generates the latex Code for one table.

    :param header: list of header cells
    :type header: listr(string)
    :param data: list of rows
    :type data: list(list(str))
    :return: The Latex code of the table
    :rtype: string
    """
    colcount = len(header)
    s = "\\begin{tabular}{"
    for i in range(colcount):
        s = s + "|c"
    s = s + "|}\n"

    # header
    s = s + "\\hline\n" + " & ".join(header) + "\\\\ \n" + "\\hline\n"

    for row in data:
        s = s + " & ".join(row) + "\\\\ \n"
        # s = s + "$" + "$ & $".join(row) + "$\\\\ \n"

    s = s + "\\hline\n\\end{tabular}\n\n"
    return s


def wrhairy_vs_dim_poly(l, h, w):
    max_vertices = 2*l-2+h

def vs_dim_formatted(vs):
    if not vs.is_valid():
        return "-"
    if not vs.exists_basis_file():
        return "?"
    return str(vs.get_dimension())

def ops_formatted(op):
    if not op.is_valid():
        return "-"
    if not op.exists_matrix_file():
        return "?"
    if not op.exists_rank_file():
        return "ok R?"
    # check if rank is mod or Q
    rank_dict = op._load_rank_dict
    if "sage_integer" in rank_dict
    return "ok R {str(op.get_matrix_rank())}"


def create_wrhairy_vs_table(v_range, l_range, h_range, w_range):
    s = ""

    header = ["l,v"] + [ str(v) for v in v_range ]
    
    for w in w_range:
        for h in h_range:
            s = s + f"\n\\smallskip\n{w} omegas, {h} hairs \n\n"
            data = []
            for l in l_range:
                data.append( [ str(l) ] + [vs_dim_formatted(WRHairyGraphComplex.WRHairyGraphVS(v,l,h,w)) for v in v_range] ) 
            s = s+latex_table(header, data)
    return s



def create_wrhairy_ops_table(v_range, l_range, h_range, w_range):
    s = ""

    header = ["l,v"] + [ str(v) for v in v_range ]
    
    for w in w_range:
        for h in h_range:
            s = s + f"\n\\smallskip\n{w} omegas, {h} hairs \n\n"
            data = []
            for l in l_range:
                data.append( [ str(l) ] + [vs_dim_formatted(WRHairyGraphComplex.WRHairyGraphVS(v,l,h,w)) for v in v_range] ) 
            s = s+latex_table(header, data)
    return s

def write_tables():
    # Generate tables
    s = create_wrhairy_vs_table(range(25),range(9), range(6), range(1,3))
    with open(latexfile_wrhairy_vs, 'w') as f:
        f.write( s )

    s = create_wrhairy_ops_table(range(25),range(9), range(6), range(1,3))
    with open(latexfile_wrhairy_ops, 'w') as f:
        f.write( s )

    s = create_wrhairy_cohom_table(range(25),range(9), range(6), range(1,3))
    with open(latexfile_wrhairy_cohom, 'w') as f:
        f.write( s )


write_tables()

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
