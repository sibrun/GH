"""General parameters."""

import sys
import os

# ---- Genearal Parameters ----
# ---- Directory Names ----
data_home_dir = "gh_data/"
data_dir = data_home_dir + "data"
plots_dir = data_home_dir + "plots"
ref_data_dir = data_home_dir + "data_ref"
log_dir = data_home_dir + "log"
temp_folder = "temp"
web_dir = "web"

# ---- Epsilon Parameters ----
square_zero_test_eps = 1e-6
commute_test_eps = 1e-6

# ---- Rank Computation ----
prime = 32189   # Prime number to be used in rank computations.
# Use sage to determine the matrix rank over the integers or over a finite field
sage_rank_options = {'integer', 'mod'}
# (modulo a prime number).

# ---- Display Parameters ----
# x width of the unit squares in the cohomology dimension plots.
x_width = 0.4
# y width of the unit squares in the cohomology dimension plots.
y_width = 0.4
# Symbol to indicate that the dimension of the cohomology is zero, but the dimension of the
zero_symbol = '*'
# vector space is not zero.
# Symbol to indicate that the dimension of the vector space is zero and hence also the dimension
zero_v_symbol = '-'
# of the cohomology.

# ---- Max Values for Sorting ----
max_sort_value = sys.maxsize     # Return value if dimension or shape is unknown

# ---- Data Tracker Parameters ----
data_tracker_timeout = 0.2

# ---- Secon Layor Info ----
# Option to plot information about the sub vector spaces of a degree slice in a graded vector space.
second_info = False
