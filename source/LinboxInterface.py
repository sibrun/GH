import os
import StoreLoad
import Parameters


def rank(linbox_option, matrix_file, prime=Parameters.prime):
    if not (linbox_option in Parameters.linbox_options):
        raise ValueError('Possible options for linbox: ' + str(Parameters.linbox_options))
    rank_suffix = 'rank.txt'
    temp_rank_file = matrix_file + rank_suffix
    linbox_path = os.path.join(os.path.curdir, "source", "rank")
    if linbox_option is "rational":
        linbox_command = "%s %s %s" % (linbox_path, matrix_file, temp_rank_file)
    else:
        linbox_command = "%s %s %s %d" % (linbox_path, matrix_file, temp_rank_file, prime)
    os.system(linbox_command)
    rank = int(StoreLoad.load_line(temp_rank_file))
    StoreLoad.delete_file_and_empty_dir(temp_rank_file)
    return rank