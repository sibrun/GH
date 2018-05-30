import os
import tempfile
import Parameters


linbox_options = {"rational", "mod"}    # Linbox options for rank computation: 'rational' (exact rank over the rational
                                        # numbers), 'mod' (rank over a finite field, i.e. all calculations modulo a
                                        # prime number).


def rank(linbox_option, matrix_file, prime=Parameters.prime):
    if not (linbox_option in linbox_options):
        raise ValueError('Possible options for linbox: ' + str(linbox_options))
    linbox_path = os.path.join(os.path.curdir, "source", "rank")
    with tempfile.NamedTemporaryFile() as temp_rank_file:
        if linbox_option is "rational":
            linbox_command = "%s %s %s" % (linbox_path, matrix_file, temp_rank_file.name)
        else:
            linbox_command = "%s %s %s %d" % (linbox_path, matrix_file, temp_rank_file.name, prime)
        os.system(linbox_command)
        rank = int(temp_rank_file.read())
    return rank
