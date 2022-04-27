"""Interface to the linbox library for rank computations"""

import os
import tempfile
import Parameters


# Linbox options for rank computation: 'rational' (exact rank over the rational
linbox_options = {"rational", "mod"}
# numbers), 'mod' (rank over a finite field, i.e. all calculations modulo a
# prime number).


def rank(linbox_option, matrix_file, prime=Parameters.prime):
    """Call the linbox library to compute the rank of a matrix.

    :param linbox_option: Options for rank computation: 'rational' (exact rank over the rational numbers),
                          'mod' (rank over a finite field, i.e. all calculations modulo a prime number).
    :type linbox_option: str
    :param matrix_file: Path to the matrix file. The matrix needs to be stored in the SMS format.
    :type matrix_file: path
    :param prime: Prime number for modulo rank calculations (Default: Parameters.prime)
    :type prime: int
    :return: Matrix rank calculated by the linbox library.
    :rtype: int

    .. seealso:: - http://www.linalg.org/
                 - https://github.com/linbox-team/linbox/blob/master/examples/rank.C
                 - http://ljk.imag.fr/membres/Jean-Guillaume.Dumas/simc.html
    """
    if not (linbox_option in linbox_options):
        raise ValueError('Possible options for linbox: ' + str(linbox_options))
    linbox_path = os.path.join(os.path.curdir, "rank_exe", "rank")
    with tempfile.NamedTemporaryFile() as temp_rank_file:
        if linbox_option is "rational":
            print("Using Linbox over rationals....")
            linbox_command = "%s %s %s" % (
                linbox_path, matrix_file, temp_rank_file.name)
        elif linbox_option is "mod":
            print(f"Using Linbox mod prime {prime}....")
            linbox_command = "%s %s %s %d" % (
                linbox_path, matrix_file, temp_rank_file.name, prime)
        else:
            raise ValueError(f"Unsupported Linbox option {linbox_option}.")
        os.system(linbox_command)
        rank = int(temp_rank_file.read())
    return rank
