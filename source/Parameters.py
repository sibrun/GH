"""General parameters."""

import sys
import os

#---- Genearal Parameters ----
#---- Directory Names ----
data_dir = "data"
plots_dir = "plots"
ref_data_dir = "data_ref"
log_dir = "log"
temp_folder = "temp"

#---- Epsilon Parameters ----
square_zero_test_eps = 1e-6
commute_test_eps = 1e-6

#---- Rank Computation ----
prime = 32189   # Prime number to be used in rank computations.
sage_rank_options = {'integer', 'mod'}  # Use sage to determine the matrix rank over the integers or over a finite field
                                        # (modulo a prime number).

#---- Display Parameters ----
x_width = 0.4           # x width of the unit squares in the cohomology dimension plots.
y_width = 0.4           # y width of the unit squares in the cohomology dimension plots.
zero_symbol = '*'       # Symbol to indicate that the dimension of the cohomology is zero, but the dimension of the
                        # vector space is not zero.
zero_v_symbol = '-'     # Symbol to indicate that the dimension of the vector space is zero and hence also the dimension
                        # of the cohomology.

#---- Max Values for Sorting ----
max_sort_value = sys.maxsize     # Return value if dimension or shape is unknown

#---- Data Tracker Parameters ----
data_tracker_timeout = 0.2

#---- Secon Layor Info ----
second_info = False     # Option to plot information about the sub vector spaces of a degree slice in a graded vector space.


