import OrdinaryGraphComplex
from sage.all import *
import GraphVectorSpace
import StoreLoad
import os
from tqdm import tqdm

maskfile_ending = ".34mask"


class Valence34Mask():
    def __init__(self, vs) -> None:
        """ """
        self.vs : GraphVectorSpace.GraphVectorSpace = vs

    def maskfile(self):
        return self.vs.get_basis_file_path() + maskfile_ending

    def exists_maskfile(self):
        return os.path.exists(self.maskfile())

    def load_mask(self):
        if not self.vs.is_valid():
            return []
        strmask = StoreLoad.load_string_list(self.maskfile())
        return [int(s) for s in strmask]

    def is_graph_34valent(self, G):
        """Whether the graph has only 3- and 4-valent vertices"""
        return all( len(G[v])<=4 for v in G.vertices() )

    def _compute_mask(self):
        vsdim = self.vs.get_dimension()
        return [(1 if self.is_graph_34valent(G) else 0)
            for G in tqdm(self.vs.get_basis(), total=vsdim, desc="Computing mask...")]

    def compute_mask(self):
        if self.vs.is_valid() and self.vs.exists_basis_file() and not self.exists_maskfile():
            print("Generating mask: ", str(self.vs))
            msk = self._compute_mask()
            strmask = [str(m) for m in msk]
            StoreLoad.store_string_list(strmask, self.maskfile())

    def index_dict(self):
        """dict from all indices to bridgeless indices"""
        msk = self.load_mask()
        mski = enumerate(msk)
        ifiltered = (i for (i, m) in mski if m == 1)
        return {i: j for (j, i) in enumerate(ifiltered)}

    def filtered_dim(self):
        return sum(self.load_mask())

    def get_34dimension(self):
        if not self.vs.is_valid():
            return 0
        msk = self.load_mask()
        return sum(msk)

    def get_P34(self):
        """matrix from get_34dimension()-dimensional space to full GVS"""
        msk = self.load_mask()
        dimsrc = sum(msk)
        dimtgt = self.vs.get_dimension()
        M = matrix(dimtgt, dimsrc, sparse = True)
        j = 0
        for (i, v) in enumerate(msk):
            if v == 1:
                M[i, j] = 1
                j=j+1
        return M

    def get_P5(self):
        """matrix from full vs to 5-valent subspace"""
        msk = self.load_mask()
        dimtgt = sum(1-v for v in msk)
        dimsrc = self.vs.get_dimension() if self.vs.is_valid() else 0
        M = matrix(dimtgt, dimsrc, sparse = True)
        j = 0
        for (i, v) in enumerate(msk):
            if v == 0:
                M[j, i] = 1
                j=j+1
        return M

    def get_34index_list(self):
        """List of indices from 34 space into full space"""
        msk = self.load_mask()
        return [i for (i,v) in enumerate(msk) if v==1]

    def get_5index_list(self):
        """List of indices from 5 space into full space"""
        msk = self.load_mask()
        return [i for (i,v) in enumerate(msk) if v==0]


def get_34cohom_dim(v,l, even_e):
    """ Compute cohomology dimension ..."""
    op1 = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(v,l, even_e)
    op2 = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(v+1,l, even_e)
    fullvs = op1.domain
    fullvs2 = op2.domain

    vs34 = Valence34Mask(fullvs)
    vs342 = Valence34Mask(fullvs2)

    D34rank = 0
    if op1.is_valid():
        D = op1.get_matrix()
        # P34 = vs34.get_P34()
        # D34 = D * P34
        i34 = vs34.get_34index_list()
        D34 = D[:, i34]
        D34rank = D34.rank()

    DD34rank = 0
    DD5rank = 0
    if op2.is_valid():
        DD = op2.get_matrix()
        # PP34 = vs342.get_P34()
        ii34 = vs342.get_34index_list()
        # DD34 = DD * PP34
        DD34 = DD[:,ii34]
        DD34rank = DD34.rank()

        # P5 = vs34.get_P5()
        # DD5 = P5 * DD * PP34
        i5 = vs34.get_5index_list()
        DD5 = DD[i5, ii34]
        DD5rank = DD5.rank()


    return vs34.get_34dimension() - D34rank -DD34rank + DD5rank


# test it
even_e = True
OGC = OrdinaryGraphComplex.OrdinaryGC(range(20), range(10), even_e, ['contract'], shift_loops_minus_vertices=100)
OGC.build_basis()
OGC.build_matrix()


# # build masks
for vs in OGC.sum_vector_space.vs_list:
    Valence34Mask(vs).compute_mask()

for l in range(3,10):
    for v in range(4,2*l-1):
        vs = OrdinaryGraphComplex.OrdinaryGVS(v,l,even_e)
        print(vs.get_dimension(), Valence34Mask(vs).get_34dimension())

# for l in range(3,10):
#     for v in range(4,2*l-1):
#         print(f"l={l},v={v}, cohomdim={ get_34cohom_dim(v,l, even_e) }")


