import os
import StoreLoad


def rank(matrix_file, prime=None):
    rank_suffix = 'rank.txt'
    temp_rank_file = matrix_file + rank_suffix
    if prime is None:
        linbox_command = "./rank %s %s" % (matrix_file, temp_rank_file)
    else:
        linbox_command = "./rank %s %s %d" % (matrix_file, temp_rank_file, prime)
    os.system(linbox_command)
    rank = int(StoreLoad.load_line(temp_rank_file))
    StoreLoad.delete_file_and_empty_dir(temp_rank_file)
    return rank