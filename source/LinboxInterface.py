"""Interface to the linbox library for rank computations"""

import os
import tempfile
import Parameters
import MatrixMethods


# Linbox options for rank computation: 'rational' (exact rank over the rational
linbox_options = {"rational", "mod", "modprecond"}
# numbers), 'mod' (rank over a finite field, i.e. all calculations modulo a
# prime number), 'modpre' (same as mod, but with preconditioning the matrix).

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
    if linbox_option not in linbox_options:
        raise ValueError('Possible options for linbox: ' + str(linbox_options))
    linbox_path = os.path.join(os.path.curdir, "rank_exe", "rank")

    rankbias = 0
    with tempfile.NamedTemporaryFile() as temp_rank_file:
        if linbox_option == "rational":
            print("Using Linbox over rationals....")
            linbox_command = "%s %s %s" % (
                linbox_path, matrix_file, temp_rank_file.name)
        elif linbox_option == "mod":
            print(f"Using Linbox mod prime {prime}....")
            linbox_command = "%s %s %s %d" % (
                linbox_path, matrix_file, temp_rank_file.name, prime)
        elif linbox_option == "modprecond":
            print(f"Using Linbox mod prime {prime} with preconditioning....")
            matrix_file_precond, rankbias = MatrixMethods.precondition_file(matrix_file, ensure_m_greater_n=True)
            linbox_command = "%s %s %s %d" % (
                linbox_path, matrix_file_precond, temp_rank_file.name, prime)
        else:
            raise ValueError(f"Unsupported Linbox option {linbox_option}.")
        ret = os.system(linbox_command)
        if ret != 0:
            raise RuntimeError("Linbox rank returned a nonzero exit code.")
        rank = int(temp_rank_file.read()) + rankbias
    return rank
