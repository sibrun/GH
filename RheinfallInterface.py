import os
import StoreLoad


def rank(command, matrix_file):
    rank_suffix = 'rank.txt'
    temp_rank_file = os.path.join(matrix_file, rank_suffix)
    command = "./%s %s %s" % (command, matrix_file, temp_rank_file)
    os.system(command)
    rank = int(StoreLoad.load_line(temp_rank_file))
    StoreLoad.delete_file_and_empty_dir(temp_rank_file)
    return rank

