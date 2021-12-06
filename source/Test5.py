import unittest
import itertools
import logging
import Log
import TestGraphComplex
import OrdinaryGraphComplex
from sage.all import *


log_file = "Forested_Unittest.log"

if __name__ == '__main__':
    tt = OrdinaryGraphComplex.OrdinaryGC(
        range(0, 15), range(0, 8), False, {'contract'})
    tt.build_basis(ignore_existing_files=True, n_jobs=1, progress_bar=True)
    tt.build_matrix(ignore_existing_files=True, n_jobs=1, progress_bar=True)
