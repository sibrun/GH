import os
import StoreLoad
import Parameters


def rank(rheinfall_option, matrix_file):
    if not (rheinfall_option in Parameters.rheinfall_options):
        raise ValueError('Possible options for rheinfall: ' + str(Parameters.rheinfall_options))
    rank_suffix = 'rank.txt'
    temp_rank_file = os.path.join(matrix_file, rank_suffix)
    rheinfall_option = "./%s %s %s" % (rheinfall_option, matrix_file, temp_rank_file)
    os.system(rheinfall_option)
    rank = int(StoreLoad.load_line(temp_rank_file))
    StoreLoad.delete_file_and_empty_dir(temp_rank_file)
    return rank

