"""Interface to the rheinfall library for rank computations"""
import os
import tempfile


rheinfall_options = {"int64", "mpq", "mpz"}  # Rheinfall options for rank computation: 'int64' (rank calculation
                                             # modulo 64 bit integer), 'mpz' (exact rank over Z), 'mpq' (exact rank over
                                             # the rational numbers).


def rank(rheinfall_option, matrix_file):
    """Call the rheinfall library to compute the rank of a matrix.

    :param rheinfall_option: Options for rank computation: 'int64' (rank calculation modulo 64 bit integer),
                             'mpz' (exact rank over Z), 'mpq' (exact rank over the rational numbers).
    :type rheinfall_option: str
    :param matrix_file: Path to the matrix file. The matrix needs to be stored in the SMS format.
    :type matrix_file: path
    :return: Matrix rank calculated by the rheinfall library.
    :rtype: int

    .. seealso:: - https://github.com/riccardomurri/rheinfall/blob/master/src.c%2B%2B/examples/rank.cpp
                 - http://ljk.imag.fr/membres/Jean-Guillaume.Dumas/simc.html
    """
    if not (rheinfall_option in rheinfall_options):
        raise ValueError('Possible options for rheinfall: ' + str(rheinfall_options))
    rheinfall_path = os.path.join(os.path.curdir, "linbox_rheinfall_rank", "rank")
    with tempfile.NamedTemporaryFile() as temp_rank_file:
        rheinfall_command = "%s-%s %s %s" % (rheinfall_path, rheinfall_option, matrix_file, temp_rank_file.name)
        os.system(rheinfall_command)
        rank = int(temp_rank_file.read())
    return rank
