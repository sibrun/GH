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

latexdir = os.path.join(Parameters.plots_dir, "latex")
latexfile_wrhairy = os.path.join(latexdir, "wrhairy.tex")

StoreLoad.generate_path(latexdir)


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
    s = "\\begin{tabular}["
    for i in range(colcount):
        s = s + "c"
    s = s + "]\n"

    # header
    s = "&".join(header) + "\\\\ \n"

    for row in data:
        s = "$" + "$ & $".join(row) + "$\\\\ \n"

    s = s + "\\end{tabular}"
    return s


def wrhairy_vs_dim_poly(l, h, w):
    max_vertices = 


def create_wrhairy_vs_table(l_range, h_range, w):
    s = ""

    header = [""] + [ str(h) for h in h_range ]
    data = []
    for l in l_range:
        data.append( [ str(l) ] + [wrhairy_vs_dim_poly(l,h, w) for h in h_range] ) 



def write_tables():
    # Generate tables
    s = create_wrhairy_vs_table(range(15), range)
    with open(latexfile_wrhairy, 'w') as f:
        f.write( s )



write_tables()

 def print_dim_and_eulerchar(self):
        for w in self.w_range:
            for h in self.h_range:
                for l in self.l_range:
                    ds = [WRHairyGraphVS(v, l, h, w).get_dimension()
                          for v in self.v_range]
                    eul = sum([(1 if j % 2 == 0 else -1) *
                              d for j, d in enumerate(ds)])
                    print("Dimensions (w,h,l) ", w,
                          h, l, ":", ds, "Euler", eul)

    def print_cohomology_dim(self):
        for w in self.w_range:
            for h in self.h_range:
                for l in self.l_range:
                    cohomdict = {}
                    for v in self.v_range:
                        D1 = ContractEdgesGO.generate_operator(v, l, h, w)
                        D2 = ContractEdgesGO.generate_operator(v+1, l, h, w)
                        try:
                            d = WRHairyGraphVS(v, l, h, w).get_dimension()
                            r1 = D1.get_matrix_rank()
                            r2 = D2.get_matrix_rank()
                            cohomdict[v] = d-r1-r2
                        except:
                            pass

                    print("Cohomology Dimensions (w,h,l) ",
                          w, h, l, ":", cohomdict)
