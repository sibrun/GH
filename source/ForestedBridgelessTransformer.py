import ForestedGraphComplex
import StoreLoad
import os
import GraphOperator
import GraphVectorSpace
from tqdm import tqdm

maskfile_ending = ".blmsk"  # the masks are stored in files with this ending

if ForestedGraphComplex.use_bridgeless:
    raise RuntimeError(
        "ForestedBridgelessTransformer: Use this module only if use_bridgeless is False!")


def get_blfilepath(path):
    return path.replace("/forested/", "/forestedbl/")


"""Contains routines for transcribing existing forested data files 
(with bridges) into bridgeless form."""


class BridgelessMask():
    def __init__(self, vs) -> None:
        """ vs must be a ForestedGraphComplex or PreForestedGVS.
        WITH BRIDGES (...that is use_bridgeless must be off)"""
        self.vs : GraphVectorSpace.GraphVectorSpace = vs
        if isinstance(vs, ForestedGraphComplex.PreForestedGVS):
            self.prevs = vs
        else:
            self.prevs = vs.preVS

    def maskfile(self):
        return self.vs.get_basis_file_path() + maskfile_ending

    def exists_maskfile(self):
        return os.path.exists(self.maskfile())

    def load_mask(self):
        strmask = StoreLoad.load_string_list(self.maskfile())
        return [int(s) for s in strmask]

    def _compute_mask(self):
        vsdim = self.vs.get_dimension()
        return [(1 if self.prevs.is_bridgeless(G) else 0) 
            for G in tqdm(self.vs.get_basis(), total=vsdim, desc="Computing mask...")]

    def compute_mask(self):
        if self.vs.is_valid() and self.vs.exists_basis_file() and not self.exists_maskfile():
            print("Generating mask: ", str(self.vs))
            msk = self._compute_mask()
            strmask = [str(m) for m in msk]
            StoreLoad.store_string_list(strmask, self.maskfile())

    def create_filtered_basisfile(self):
        if not (self.vs.is_valid() and self.vs.exists_basis_file()):
            return
        newbasisfile = get_blfilepath(self.vs.get_basis_file_path())
        if os.path.exists(newbasisfile):
            return

        if not self.exists_maskfile():
            self.compute_mask()

        print("Filtering basis: ",  str(self.vs))

        mask = self.load_mask()
        basisg6 = self.vs.get_basis_g6()

        filteredbasis = [s for (s, m) in zip(basisg6, mask) if m == 1]
        filteredbasis.insert(0, str(len(filteredbasis)))
        StoreLoad.store_string_list(filteredbasis, newbasisfile)

    def index_dict(self):
        """dict from all indices to bridgeless indices"""
        msk = self.load_mask()
        mski = enumerate(msk)
        ifiltered = (i for (i, m) in mski if m == 1)
        return {i: j for (j, i) in enumerate(ifiltered)}

    def filtered_dim(self):
        return sum(self.load_mask())


class BridgeLessMaskOM():
    def __init__(self, op: GraphOperator.OperatorMatrix):
        self.op = op
        self.domain = BridgelessMask(op.domain)
        self.target = BridgelessMask(op.target)

    def create_filtered_matrixfile(self):
        if not (self.op.is_valid() and self.op.exists_matrix_file()):
            return

        newfilepath = get_blfilepath(self.op.get_matrix_file_path())
        if os.path.exists(newfilepath):
            return

        print("Translating operator: ", str(self.op))
        idict1 = self.domain.index_dict()
        idict2 = self.target.index_dict()

        lst, oldshape = self.op._load_matrix_list()
        newshape = (self.domain.filtered_dim(), self.target.filtered_dim())

        matrix_list = [(idict1[i], idict2[j], v)
                       for (i, j, v) in lst if (i in idict1 and j in idict2)]

        # copied
        data_type = "M"
        (d, t) = newshape
        stringList = []
        stringList.append("%d %d %s" % (d, t, data_type))
        for (i, j, v) in matrix_list:
            stringList.append("%d %d %d" % (i + 1, j + 1, v))
        stringList.append("0 0 0")
        StoreLoad.store_string_list(stringList, newfilepath)


h_range = range(0, 9)
v_range = range(0, 25)
l_range = range(0, 10)
m_range = range(0, 25)
for h in h_range:
    for l in l_range:
        for v in v_range:
            for m in m_range:
                bm = BridgelessMask(
                    ForestedGraphComplex.PreForestedGVS(v, l, m, h))
                bm.create_filtered_basisfile()
                for even_edges in [True]: # [True, False]: # odd not supported due to tadpoles
                    bm = BridgelessMask(
                        ForestedGraphComplex.ForestedGVS(v, l, m, h, even_edges))
                    bm.create_filtered_basisfile()


for h in h_range:
    for l in l_range:
        for v in v_range:
            for m in m_range:
                for even_edges in [True]: # [True, False]: # odd not supported due to tadpoles
                    bo = BridgeLessMaskOM(
                        ForestedGraphComplex.ContractEdgesGO.generate_operator(v, l, m, h, even_edges))
                    bo.create_filtered_matrixfile()
                    bo = BridgeLessMaskOM(
                        ForestedGraphComplex.UnmarkEdgesGO.generate_operator(v, l, m, h, even_edges))
                    bo.create_filtered_matrixfile()
