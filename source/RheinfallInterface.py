import os
import tempfile
import StoreLoad
import Parameters


rheinfall_options = {"int64", "mpq", "mpz"}  # Rheinfall options for rank computation: 'int64' (rank calculation
                                             # modulo 64 bit integer), 'mpz' (exact rank over Z), 'mpq' (exact rank over
                                             # the rational numbers).

def rank(rheinfall_option, matrix_file):
    if not (rheinfall_option in rheinfall_options):
        raise ValueError('Possible options for rheinfall: ' + str(rheinfall_options))
    #rank_suffix = 'rank.txt'
    #temp_rank_file = matrix_file + rank_suffix
    rheinfall_path = os.path.join(os.path.curdir, "source", "rank")
    with tempfile.NamedTemporaryFile() as temp_rank_file:
        rheinfall_command = "%s-%s %s %s" % (rheinfall_path, rheinfall_option, matrix_file, temp_rank_file)
        os.system(rheinfall_command)
        rank = int(temp_rank_file.read())
    StoreLoad.delete_file_and_empty_dir(temp_rank_file)
    return rank

