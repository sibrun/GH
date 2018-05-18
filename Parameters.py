"""General parameters."""

import sys

#---- Genearal Parameters ----
#---- Directory Names ----
data_dir = "data"
plots_dir = "plots"
ref_data_dir = "data_ref"
log_dir = "log"

#---- Epsilon Parameters ----
square_zero_test_eps = 1e-6
commute_test_eps = 1e-6
estimate_rank_eps = 1e-6

#---- Rank Computation ----
primes = [32189, 31277, 32183, 31121]   # List of prime numbers to be used in rank computations.
min_size_for_rank_estimate = 40     # Minimal size for the number of rows or columns of the matrix to allow
                                    # estimating the rank. The method is not reliable for small matrices.

#---- Display Parameters ----
x_width = 0.4           # x width of the unit squares in the cohomology dimension plots.
y_width = 0.4           # y width of the unit squares in the cohomology dimension plots.
zero_symbol = '*'       # Symbol to indicate that the dimension of the cohomology is zero, but the dimension of the
                        # vector space is not zero.
zero_v_symbol = '-'     # Symbol to indicate that the dimension of the vector space is zero and hence also the dimension
                        # of the cohomology.

#---- Max Values for Sorting ----
max_sort_value = sys.maxint     # Return value if dimension or shape is unknown

#---- Data Tracker Parameters ----
data_tracker_timeout = 0.2

#---- Secon Layor Info ----
second_info = False     # Option to plot information about the sub vector spaces of a degree slice in a graded vector space.

#---- Global Variables ----
zero_hairs = False      # Option to include zero hairs in the hairy graph complexes.
