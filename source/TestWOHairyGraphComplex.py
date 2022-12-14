import unittest
import itertools
import logging
import Log
import TestGraphComplex
from WOHairyGraphComplex import *
from sage.all import *


# GC = WOHairyGC(range(5), range(3), range(2), range(3), ['contract'])


# GC.build_basis()

VS1 = WOHairyGraphVS(1, 3, 0, 4)
print(VS1.is_valid())
VS1.build_basis(ignore_existing_files=True)
VS1.display_basis_plots()
