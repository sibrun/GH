import argparse
import logging
import os
import Profiling
import StoreLoad as SL
import OrdinaryGraphComplex as OGC


log_dir = 'log'


def positive_int(value):
    value = int(value)
    if value <= 0:
        raise argparse.ArgumentTypeError('positive integer expected')
    return value


def range_type(arg):
    (min, max) = map(positive_int, arg.split(','))
    if min >= max:
        raise argparse.ArgumentTypeError('range min,max with 0 < min < max expected')
    return range(min, max)


graph_types = ['ogc']

log_levels = ['info', 'warn', 'error']

parser = argparse.ArgumentParser(description='Compute the homology of a graph complex')

parser.add_argument('graph_type', type=str, choices=graph_types, help='type of the graph complex')
parser.add_argument('-even_edges', action='store_true')
parser.add_argument('-odd_edges', action='store_true')
parser.add_argument('-v_range', type=range_type, help='range min,max for number of vertices')
parser.add_argument('-l_range', type=range_type, help='range min,max for number of loops')
parser.add_argument('-ignore_existing', action='store_true', help='ignore existing files')
parser.add_argument('-n_jobs', type=positive_int, default=1, help='number of parallel processes')
parser.add_argument('-profile', action='store_true', help='profiling')
parser.add_argument('-log', type=str, choices=log_levels, help='logging level')
parser.add_argument('-build', action='store_true', help='just build vector space basis and operator matrix')
parser.add_argument('-ranks', action='store_true', help='just compute matrix ranks')
parser.add_argument('-cohomology', action='store_true', help='just compute cohomology')

args = parser.parse_args()


@Profiling.cond_decorator(args.profile, Profiling.profile(log_dir))
def main(graph_complex):
    graph_complex.build(ignore_existing_files=args.ignore_existing, n_jobs=args.n_jobs)
    graph_complex.compute_ranks(ignore_existing_files=args.ignore_existing)
    graph_complex.compute_cohomology_dim()
    graph_complex.plot_cohomology_dim()


@Profiling.cond_decorator(args.profile, Profiling.profile(log_dir))
def build(graph_complex):
    graph_complex.build(ignore_existing_files=args.ignore_existing, n_jobs=args.n_jobs)


@Profiling.cond_decorator(args.profile, Profiling.profile(log_dir))
def ranks(graph_complex):
    graph_complex.compute_ranks(ignore_existing_files=args.ignore_existing)


@Profiling.cond_decorator(args.profile, Profiling.profile(log_dir))
def cohomology(graph_complex):
    graph_complex.compute_cohomology_dim()
    graph_complex.plot_cohomology_dim()


class MissingArgumentError(RuntimeError):
    pass


if __name__ == "__main__":
    if args.log is not None:
        log_file = args.graph_type + '.log'
        log_path = os.path.join(log_dir, log_file)
        SL.generate_path(log_path)
        logging.basicConfig(filename=log_path, level=args.log.upper())
        logging.warn("\n###########################\n" + "----- Graph Homology -----")

    graph_complex = None

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

        graph_complex = OGC.OrdinaryGC(args.v_range, args.l_range, even_edges)

    if graph_complex is None:
        raise ValueError('graph complex not specified')

    if not (args.build or args.ranks or args.cohomology):
        main(graph_complex)
    else:
        if args.build:
            build(graph_complex)
        if args.ranks:
            ranks(graph_complex)
        if args.cohomology:
            cohomology(graph_complex)