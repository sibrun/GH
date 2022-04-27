import unittest
import itertools
import logging
import Log
import TestGraphComplex
import ForestedGraphComplex
from sage.all import *

VS = ForestedGraphComplex.ForestedGVS(6, 4, 3, 0, True)
VS.build_basis()
