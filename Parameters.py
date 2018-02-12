import sys

#---- Genearal Parameters ----
#---- Directory Names ----
data_dir = "data"
plots_dir = "plots"
ref_data_dir = "data_ref"
log_dir = "log"
temp_folder = 'temp'

#---- Test Parameters ----
eps = 1e-6

#---- Display Parameters ----
x_width = 0.7
y_width = 0.7
zero_symbol = '*'
x_plots = 2

#---- Progress Bar Parameters ----
pbar_steps = 200
timeout = 0.1
# Minimum progress display update interval, in iterations.
# If 0 and ``dynamic_miniters``, will automatically adjust to equal
# ``mininterval`` (more CPU efficient, good for tight loops).
# If > 0, will skip display of specified number of iterations.


#---- Max Values for Sorting ----
MAX_ENTRIES = sys.maxint  # return value if number of entries is unkown, i.e. if no matrix file
MAX_DIMENSION = sys.maxint  # return value if dimension is unkown, i.e. if no basis file
