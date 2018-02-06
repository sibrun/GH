import argparse
import logging
import os
import Profiling
import StoreLoad as SL
import OrdinaryGraphComplex as OGC
import HairyGraphComplex as HGC
import Parameters


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


graph_types = ['ogc', 'hgc']

log_levels = ['info', 'warn', 'error']

parser = argparse.ArgumentParser(description='Compute the homology of a graph complex')

parser.add_argument('graph_type', type=str, choices=graph_types, help='type of the graph complex')
parser.add_argument('-even_e', action='store_true', help='even edges')
parser.add_argument('-odd_e', action='store_true', help='odd edges')
parser.add_argument('-even_h', action='store_true', help='even hairs')
parser.add_argument('-odd_h', action='store_true', help='odd edges')
parser.add_argument('-v', type=range_type, help='range min,max for number of vertices')
parser.add_argument('-l', type=range_type, help='range min,max for number of loops')
parser.add_argument('-hairs', type=range_type, help='range min,max for number of hairs')
parser.add_argument('-ignore_ex', action='store_true', help='ignore existing files')
parser.add_argument('-n_jobs', type=positive_int, default=1, help='number of parallel processes')
parser.add_argument('-profile', action='store_true', help='profiling')
parser.add_argument('-log', type=str, choices=log_levels, help='logging level')
parser.add_argument('-build', action='store_true', help='just build vector space basis and operator matrix')
parser.add_argument('-rank', action='store_true', help='just compute matrix ranks')
parser.add_argument('-coho', action='store_true', help='just compute cohomology')

args = parser.parse_args()


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def main(graph_complex):
    graph_complex.build(ignore_existing_files=args.ignore_ex, n_jobs=args.n_jobs)
    graph_complex.sort_member(work_estimate=False)
    graph_complex.compute_ranks(ignore_existing_files=args.ignore_ex)
    graph_complex.get_cohomology_dim()
    graph_complex.plot_cohomology_dim()


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def build(graph_complex):
    graph_complex.build(ignore_existing_files=args.ignore_ex, n_jobs=args.n_jobs)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def rank(graph_complex):
    graph_complex.sort_member(work_estimate=False)
    graph_complex.compute_ranks(ignore_existing_files=args.ignore_ex)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def cohomology(graph_complex):
    graph_complex.sort_member(work_estimate=False)
    graph_complex.get_cohomology_dim()
    graph_complex.plot_cohomology_dim()


class MissingArgumentError(RuntimeError):
    pass


if __name__ == "__main__":
    if args.log is not None:
        log_file = args.graph_type + '.log'
        log_path = os.path.join(Parameters.log_dir, log_file)
        SL.generate_path(log_path)
        logging.basicConfig(filename=log_path, level=args.log.upper())
        logging.warn("\n###########################\n" + "----- Graph Homology -----")

    graph_complex = None

    if args.graph_type in {'ogc', 'hgc'}:
        if args.even_e:
                even_edges = True
        elif args.odd_e:
                even_edges = False
        else:
            raise MissingArgumentError('specify -even_e or -odd_e')
        if args.v is None:
            raise MissingArgumentError('specify -v: range for number of vertices')
        if args.l is None:
            raise MissingArgumentError('specify -l: range for number of loops')

    if args.graph_type == 'hgc':
        if args.even_h:
                even_hairs = True
        elif args.odd_e:
                even_hairs = False
        else:
            raise MissingArgumentError('specify -even_h or -odd_h')
        if args.hairs is None:
            raise MissingArgumentError('specify -hairs: range for number of hairs')

    if args.graph_type == 'ogc':
        graph_complex = OGC.OrdinaryGC(args.v, args.l, even_edges)

    elif args.graph_type == 'hgc':
        graph_complex = HGC.HairyGC(args.v, args.l, args.hairs, even_edges, even_hairs)

    if graph_complex is None:
        raise ValueError('graph complex not specified')

    if not (args.build or args.rank or args.coho):
        main(graph_complex)
    else:
        if args.build:
            build(graph_complex)
        if args.rank:
            rank(graph_complex)
        if args.coho:
            cohomology(graph_complex)
