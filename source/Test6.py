import unittest
import itertools
import logging
import Log
import TestGraphComplex
import OrdinaryGraphComplex
from sage.all import *
import StoreLoad
import OrdinaryMerkulovComplex

even_e = True
OGC = OrdinaryMerkulovComplex.OrdinaryMerkulovGC(range(20), range(9), even_e, ['contract'])
OGC.build_basis()
OGC.build_matrix()
OGC.compute_rank(sage="mod")

OGC.print_cohom()