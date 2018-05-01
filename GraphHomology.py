import argparse
import multiprocessing as mp
import Log
import Profiling
import OrdinaryGraphComplex as OGC
import HairyGraphComplex as HGC
import Parameters


logger = Log.logger.getChild('main')


def positive_int(value):
    value = int(value)
    if value <= 0:
        raise argparse.ArgumentTypeError('positive integer expected')
    return value

def non_negative_int(value):
    value = int(value)
    if value < 0:
        raise argparse.ArgumentTypeError('non negative integer expected')
    return value


def non_negative_range_type(arg):
    (min, max) = map(non_negative_int, arg.split(','))
    if min >= max:
        raise argparse.ArgumentTypeError('range min,max with 0 < min < max expected')
    return range(min, max)

graph_types = ['ordinary', 'hairy']
differentials = ['contract', 'et1h']

parser = argparse.ArgumentParser(description='Compute the homology of a graph complex')

parser.add_argument('graph_type', type=str, choices=graph_types, help='type of the graphs')
parser.add_argument('dif1', type=str, choices=differentials, help='differential 1')
parser.add_argument('-dif2', type=str, choices=differentials, default=None, help='differential 2')
parser.add_argument('-even_e', action='store_true', help='even edges')
parser.add_argument('-odd_e', action='store_true', help='odd edges')
parser.add_argument('-even_h', action='store_true', help='even hairs')
parser.add_argument('-odd_h', action='store_true', help='odd edges')
parser.add_argument('-v', type=non_negative_range_type, help='range min,max for number of vertices')
parser.add_argument('-l', type=non_negative_range_type, help='range min,max for number of loops')
parser.add_argument('-hairs', type=non_negative_range_type, help='range min,max for number of hairs')
parser.add_argument('-ignore_ex', action='store_true', help='ignore existing files')
parser.add_argument('-n_jobs', type=positive_int, default=1, help='number of parallel processes')
parser.add_argument('-pbar', action='store_true', help='show progressbar')
parser.add_argument('-profile', action='store_true', help='profiling')
parser.add_argument('-log', type=str, choices=Log.log_levels_dict.keys(), help='logging level')
parser.add_argument('-info', action='store_true', help='display info during calculations in browser')
parser.add_argument('-exact', action='store_true', help='exact matrix rank computation')
parser.add_argument('-n_primes', type=non_negative_int, default=1, help='compute matrix rank modulo n_primes different prime numbers')
parser.add_argument('-est', action='store_true', help="estimate matrix rank")
parser.add_argument('-build', action='store_true', help='build vector space basis and operator matrix')
parser.add_argument('-build_b', action='store_true', help='build vector space basis')
parser.add_argument('-build_op', action='store_true', help='build operator matrix')
parser.add_argument('-rank', action='store_true', help='compute matrix ranks')
parser.add_argument('-cohomology', action='store_true', help='compute cohomology')
parser.add_argument('-square_zero_test', action='store_true', help='square zero test')
parser.add_argument('-commutativity_test', action='store_true', help='test commutativity of differentials')
parser.add_argument('-anti_commute', action='store_true', help='test commutativity of differentials')

args = parser.parse_args()


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def build_basis(graph_complex):
    graph_complex.build_basis(ignore_existing_files=args.ignore_ex, n_jobs=args.n_jobs, progress_bar=args.pbar,
                              info_tracker=args.info)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def build_operator(graph_complex):
    graph_complex.build_matrix(ignore_existing_files=args.ignore_ex, n_jobs=args.n_jobs, progress_bar=args.pbar,
                               info_tracker=args.info)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def square_zero_test(graph_complex):
    graph_complex.square_zero_test()


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def test_commutativity(graph_complex):
    graph_complex.test_pairwise_commutativity(anti_commute=args.anti_commute)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def rank(graph_complex):
    graph_complex.compute_rank(exact=args.exact, n_primes=args.n_primes, estimate=args.est,
                               ignore_existing_files=args.ignore_ex, n_jobs=args.n_jobs, info_tracker=args.info)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def cohomology(graph_complex):
    graph_complex.plot_cohomology_dim()


class MissingArgumentError(RuntimeError):
    pass


if __name__ == "__main__":
    if args.log is not None:
        complex = args.graph_type + args.dif1
        if args.dif2 is not None:
            complex += args.dif2
        Log.set_log_level(args.log)
        log_file = complex + '.log'
        Log.set_log_file(log_file)

    logger.warn("\n###########################\n" + "----- Graph Homology -----")

    graph_complex = None

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

    if args.graph_type == 'ordinary':
        if not args.dif1 == 'contract':
            raise ValueError('only contract edges differential implemented for ordinary graphs')
        graph_complex = OGC.OrdinaryGC(args.v, args.l, even_edges)

    if args.graph_type == 'hairy':
        if args.even_h:
                even_hairs = True
        elif args.odd_h:
                even_hairs = False
        else:
            raise MissingArgumentError('specify -even_h or -odd_h')

        if args.hairs is None:
            raise MissingArgumentError('specify -hairs: range for number of hairs')

        differentials = [args.dif1]
        if args.dif2 is not None:
            differentials.append(args.dif2)

        graph_complex = HGC.HairyGC(args.v, args.l, args.hairs, even_edges, even_hairs, differentials)

    if args.build_b:
        build_basis(graph_complex)
    if args.build_op:
        build_operator(graph_complex)
    if args.square_zero_test:
        square_zero_test(graph_complex)
    if args.commutativity_test:
        test_commutativity(graph_complex)
    if args.rank:
        rank(graph_complex)
    if args.cohomology:
        cohomology(graph_complex)
