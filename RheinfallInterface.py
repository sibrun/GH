import os


def rank_int64(matrix_file, rank_file):
    command = "./rank-int64 %s %s" % (matrix_file, rank_file)
    os.system(command)