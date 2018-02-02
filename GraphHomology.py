import argparse
import logging
import os
import OrdinaryGraphComplex as OGC
import Profiling
import StoreLoad as SL

reload(OGC)
reload(Profiling)
reload(SL)

log_dir = 'log'


def positive_int(value):
    if type(value) is not int or value <= 0:
        raise argparse.ArgumentTypeError('positive integer value expected')
    return value


supported_graph_types = ['ogc']

log_levels = ['info', 'warn', 'error']

parser = argparse.ArgumentParser(description='Compute the homology of a graph complex')

parser.add_argument('graph_type', type=str, choices=supported_graph_types, help='type of the graph complex')
parser.add_argument('-even_edges', action='store_true')
parser.add_argument('-odd_edges', action='store_true')
parser.add_argument('-v_range',  nargs=2, type=int, help='range for number of vertices')
parser.add_argument('-l_range', nargs=2, type=int, help='range for number of loops')
parser.add_argument('-ignore_existing', action='store_true', help='ignore existing files')
parser.add_argument('-n_jobs', type=int, default=1, help='number of parallel processes')
parser.add_argument('-profile', action='store_true')
parser.add_argument('-log', type=str, choices=log_levels, help='logging level')

args = parser.parse_args()


@Profiling.cond_decorator(args.profile, Profiling.profile(log_dir))
def ogc_main(v_range, l_range, even_edges):
    n_jobs = args.n_jobs
    ogc = OGC.OrdinaryGC(v_range, l_range, even_edges)
    ogc.build_basis(ignore_existing_files=args.ignore_existing)
    ogc.build_operator_matrix(ignore_existing_files=args.ignore_existing, n_jobs=n_jobs)
    ogc.compute_ranks(ignore_existing_files=args.ignore_existing)
    ogc.compute_cohomology_dim()
    ogc.plot_cohomology_dim()


class MissingArgumentError(RuntimeError):
    pass


if __name__ == "__main__":
    if args.log is not None:
        log_file = args.graph_type + '.log'
        log_path = os.path.join(log_dir, log_file)
        SL.generate_path(log_path)
        logging.basicConfig(filename=log_path, level=args.log.upper())
        logging.warn("\n###########################\n" + "----- Graph Homology -----")

    if args.graph_type == 'ogc':
        if args.even_edges:
                even_edges = True
        elif args.odd_edges:
                even_edges = False
        else:
            raise MissingArgumentError('specify even_edges or odd_edges')
        if args.v_range is None:
            raise MissingArgumentError('specify v_range: range for number of vertices')
        if args.l_range is None:
            raise MissingArgumentError('specify v_range: range for number of vertices')
        v_range = range(args.v_range[0], args.v_range[1])
        l_range = range(args.l_range[0], args.l_range[1])
        ogc_main(v_range, l_range, even_edges)