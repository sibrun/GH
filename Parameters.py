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
estimate_rank_eps = 1e-3

#---- Primes ----
primes = [3036995833, 3036996247, 3036996491, 3036997217, 3036997631, 3036997933]

#---- Display Parameters ----
x_width = 0.7
y_width = 0.7
zero_symbol = '*'
x_plots = 2

#---- Progress Bar Parameters ----
pbar_steps = 200
timeout = 0.1

#---- Max Values for Sorting ----
MAX_ENTRIES = sys.maxint  # return value if number of entries is unknown, i.e. if no matrix file
MAX_DIMENSION = sys.maxint  # return value if dimension is unknown, i.e. if no basis file
