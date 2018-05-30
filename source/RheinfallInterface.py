import os
import tempfile


rheinfall_options = {"int64", "mpq", "mpz"}  # Rheinfall options for rank computation: 'int64' (rank calculation
                                             # modulo 64 bit integer), 'mpz' (exact rank over Z), 'mpq' (exact rank over
                                             # the rational numbers).

def rank(rheinfall_option, matrix_file):
    if not (rheinfall_option in rheinfall_options):
        raise ValueError('Possible options for rheinfall: ' + str(rheinfall_options))
    rheinfall_path = os.path.join(os.path.curdir, "source", "rank")
    with tempfile.NamedTemporaryFile() as temp_rank_file:
        rheinfall_command = "%s-%s %s %s" % (rheinfall_path, rheinfall_option, matrix_file, temp_rank_file.name)
        os.system(rheinfall_command)
        rank = int(temp_rank_file.read())
    return rank
