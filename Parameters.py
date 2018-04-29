import sys

#---- Genearal Parameters ----
#---- Directory Names ----
data_dir = "data"
plots_dir = "plots"
ref_data_dir = "data_ref"
log_dir = "log"
temp_folder = 'temp'

#---- Eps Parameters ----
square_zero_test_eps = 1e-6
commute_test_eps = 1e-6
estimate_rank_eps = 1e-6

#---- Rank Computation ----
primes = [32189, 31277, 32183, 31121, 32143, 31481, 32233]
min_size_for_rank_estimate = 40

#---- Display Parameters ----
x_width = 0.7
y_width = 0.7
zero_symbol = '*'
x_plots = 2

#---- Progress Bar Parameters ----
pbar_steps = 200
pbar_timeout = 0.1

#---- Max Values for Sorting ----
max_sort_value = sys.maxint  # return value if dimension or shape is unknown

#---- Data Tracker Parameters ----
data_tracker_timeout = 0.2
