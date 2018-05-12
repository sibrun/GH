"""Compute the cohomology dimensions of graph complexes as well as bicomplexes.

Build graph complexes as well as bicomplexes.
Test whether an operator squares to zero, i.e. is a differential.
Test whether two operators anti-commute or commute, i.e. build a differential for a bicomplex.
Plot the cohomology dimensions of a graph complex.

This module parses command line arguments and contains the main function.
Building the basis and operator matrices as well as computing the matrix rank can be done in parallel.

There are options to ignore existing files (-ignore_ex), to display informations (-info), to show a progress bar (-pbar),
for logging (-log warning), and for profiling (-profile).

Display informations and show a progress bar are only available if not several processes work in parallel
(default: -n_jobs=1).

Prerequisits:
    sagemath: This software is based on the sage math library.
        Download sagemath from:

            http://www.sagemath.org

        Activate a sage shell with

            $ sage -sh

            or use the command $ sage --python ... instead of $ python ... .

    tqdm: Install the tqdm module in the sage python version:

            (sage-sh)$ pip -install tqdm

    pandas: Install the tqdm module in the sage python version:

            (sage-sh)$ pip -install pandas

Examples:
    Ordinary graph complex:
        Build the basis:

            $ python GraphHomology.py ordinary -op1 contract -v 3,12 -l 0,9 -odd_e -n_jobs 4 -build_b

        Build the operator matrices for the differentials 'contract edges' and 'delete edges':

            $ python GraphHomology.py ordinary -op1 contract -v 3,12 -l 0,9 -odd_e -n_jobs 4 -build_op
            $ python GraphHomology.py ordinary -op1 delete_e -v 3,12 -l 0,9 -odd_e -n_jobs 4 -build_op

        Compute the ranks of the operator matrices:

            $ python GraphHomology.py ordinary -op1 contract -v 3,12 -l 0,9 -odd_e -n_jobs 4 -n_primes 1 -rank
            $ python GraphHomology.py ordinary -op1 delete_e -v 3,12 -l 0,9 -odd_e -n_jobs 4 -n_primes 1 -rank

        Test whether the operators square to zero, i.e. build a differential, test whether the operators anti-commute, and
            plot the cohomology dimensions of the respective graph complexes:

            $ python GraphHomology.py ordinary -op1 contract -op2 delete_e -v 0,12 -l 0,9 -odd_e -square_zero -anti_commute -cohomology

    Ordinary graph bicomplex with the differentials 'contract edges' and 'delete edges':

        $ python GraphHomology.py ordinary -bicomplex ce_dele -d 0,18 -odd_e -n_jobs 4 -build_b -build_op -rank -square_zero -cohomology

    Hairy graph complex:
        Build the basis:

            $ python GraphHomology.py hairy -op1 contract -v 3,11 -l 0,10 -hairs 0,9 -odd_e -odd_h -n_jobs 4 -build_b

        Build the operator matrices for the differentials 'contract edges' and 'edge to one hair':

            $ python GraphHomology.py hairy -op1 contract -v 3,11 -l 0,10 -hairs 0,9 -odd_e -odd_h -n_jobs 4 -build_op
            $ python GraphHomology.py hairy -op1 et1h -v 3,11 -l 0,10 -hairs 0,9 -odd_e -odd_h -n_jobs 4 -build_op

        Compute the ranks of the operator matrices:

            $ python GraphHomology.py hairy -op1 contract -v 3,11 -l 0,10 -hairs 0,9 -odd_e -odd_h -n_jobs 4 -n_primes 1 -rank
            $ python GraphHomology.py hairy -op1 et1h -v 3,11 -l 0,10 -hairs 0,9 -odd_e -odd_h -n_jobs 4 -n_primes 1 -rank

        Test whether the operators square to zero, i.e. build a differential, test whether the operators anti-commute, and
            plot the cohomology dimensions of the respective graph complexes:

            $ python GraphHomology.py hairy -op1 contract -op2 et1h -v 3,11 -l 0,10 -hairs 0,9 -odd_e -odd_h -square_zero -anti_commute -cohomology

    Hairy graph bicomplex with the differentials 'contract edges' and 'edge to one hair':

        $ python GraphHomology.py hairy -bicomplex ce_et1h -d 0,14 -h_min ,-12,-1 -odd_e -odd_h -n_jobs 4 -build_b -build_op -rank -square_zero -cohomology
"""

import argparse
import Log
import Profiling
import OrdinaryGraphComplex
import OrdinaryGraphBiComplex
import HairyGraphComplex
import HairyGraphBiComplex
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


def range_type(arg):
    (prefix, min, max) = arg.split(',')
    (min, max) = map(int, (min, max))
    if min >= max:
        raise argparse.ArgumentTypeError('range min,max with min < max expected')
    return range(min, max)


graph_types = ['ordinary', 'hairy']
operators = ['contract', 'delete_e', 'et1h']
bicomplexes = ['ce_et1h', 'ce_dele']

parser = argparse.ArgumentParser(description='Compute the homology of a graph complex')

parser.add_argument('graph_type', type=str, choices=graph_types, help='type of the graphs')
parser.add_argument('-bicomplex', type=str, choices=bicomplexes, help='bicomplex')
parser.add_argument('-op1', type=str, choices=operators, help='operator 1')
parser.add_argument('-op2', type=str, choices=operators, default=None, help='operator 2')
parser.add_argument('-even_e', action='store_true', help='even edges')
parser.add_argument('-odd_e', action='store_true', help='odd edges')
parser.add_argument('-even_h', action='store_true', help='even hairs')
parser.add_argument('-odd_h', action='store_true', help='odd hairs')
parser.add_argument('-v', type=non_negative_range_type, help='range min,max for number of vertices')
parser.add_argument('-l', type=non_negative_range_type, help='range min,max for number of loops')
parser.add_argument('-hairs', type=non_negative_range_type, help='range min,max for number of hairs')
parser.add_argument('-d', type=non_negative_range_type, help='range min,max for degree of degree slices in bicomplex')
parser.add_argument('-h_min', type=range_type, help='range min,max for minimal number of hairs in a degree slice of a bicomplex')
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
parser.add_argument('-square_zero', action='store_true', help='square zero test')
parser.add_argument('-anti_commute', action='store_true', help='test anti-commutativity of differentials')
parser.add_argument('-commute', action='store_true', help='test commutativity of differentials')

args = parser.parse_args()


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def build_basis(graph_complex):
    logger.warn("\n----- Build Basis -----\n")
    graph_complex.build_basis(ignore_existing_files=args.ignore_ex, n_jobs=args.n_jobs, progress_bar=args.pbar,
                              info_tracker=args.info)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def build_operator(graph_complex):
    logger.warn("\n----- Build Matrix -----\n")
    graph_complex.build_matrix(ignore_existing_files=args.ignore_ex, n_jobs=args.n_jobs, progress_bar=args.pbar,
                               info_tracker=args.info)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def square_zero_test(graph_complex):
    logger.warn("\n----- Square Zero Test -----\n")
    graph_complex.square_zero_test()


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def test_anti_commutativity(graph_complex):
    logger.warn("\n----- Anti-commutativity Test -----\n")
    graph_complex.test_pairwise_anti_commutativity(commute=False)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def test_commutativity(graph_complex):
    logger.warn("\n----- Commutativity Test -----\n")
    graph_complex.test_pairwise_anti_commutativity(commute=True)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def rank(graph_complex):
    logger.warn("\n----- Compute Ranks -----\n")
    graph_complex.compute_rank(exact=args.exact, n_primes=args.n_primes, estimate=args.est,
                               ignore_existing_files=args.ignore_ex, n_jobs=args.n_jobs, info_tracker=args.info)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def cohomology(graph_complex):
    logger.warn("\n----- Compute Cohomology -----\n")
    graph_complex.plot_cohomology_dim()


class MissingArgumentError(RuntimeError):
    pass


if __name__ == "__main__":
    if args.log is not None:
        complex = args.graph_type
        if args.bicomplex is not None:
            complex += '_' + args.bicomplex
        if args.op1 is not None:
            complex += '_' + args.op1
        if args.op2 is not None:
            complex += '_' + args.op2
        Log.set_log_level(args.log)
        log_file = complex + '.log'
        Log.set_log_file(log_file)

    logger.warn("\n###########################\n" + "----- Graph Homology -----")

    operators = []
    if args.op1 is not None:
        operators.append(args.op1)
    if args.op2 is not None:
        operators.append(args.op2)

    if args.even_e:
            even_edges = True
    elif args.odd_e:
            even_edges = False
    else:
        raise MissingArgumentError('specify -even_e or -odd_e')

    if args.graph_type == 'ordinary':
        if args.bicomplex is not None:
            if args.bicomplex == 'ce_dele':
                if args.d is None:
                    raise MissingArgumentError('specify -d: range for degree of degree slices in bicomplex')

                graph_complex = OrdinaryGraphBiComplex.OrdinaryCeDeleBiGC(args.d, args.even_e)
            else:
                raise ValueError('Ordinary graphs bicomplexes: ce_dele')

        elif len(operators) > 0 and set(operators) <= {'contract', 'delete_e'}:
            if args.v is None:
                raise MissingArgumentError('specify -v: range for number of vertices')
            if args.l is None:
                raise MissingArgumentError('specify -l: range for number of loops')

            graph_complex = OrdinaryGraphComplex.OrdinaryGC(args.v, args.l, even_edges, operators)
        else:
            raise ValueError('Differentials for ordinary graph complex: contract, delete_e')

    elif args.graph_type == 'hairy':
        if args.even_h:
                even_hairs = True
        elif args.odd_h:
                even_hairs = False
        else:
            raise MissingArgumentError('specify -even_h or -odd_h')

        if args.bicomplex is not None:
            if args.bicomplex == 'ce_et1h':
                if args.d is None:
                    raise MissingArgumentError('specify -d: range for degree of degree slices in bicomplex')
                if args.h_min is None:
                    raise MissingArgumentError('specify -h_min: range for minimal number of hairs of degree slices in bicomplex')

                graph_complex = HairyGraphBiComplex.HairyCeEt1hBiGC(args.d, args.h_min, args.even_e, args.even_h)
            else:
                raise ValueError('Hairy graphs bicomplexes: ce_et1h')

        elif len(operators) > 0 and set(operators) <= {'contract', 'et1h'}:
            if args.v is None:
                raise MissingArgumentError('specify -v: range for number of vertices')
            if args.l is None:
                raise MissingArgumentError('specify -l: range for number of loops')
            if args.hairs is None:
                raise MissingArgumentError('specify -hairs: range for number of hairs')

            graph_complex = HairyGraphComplex.HairyGC(args.v, args.l, args.hairs, even_edges, even_hairs, operators)
        else:
            raise ValueError('Differentials for hairy graph complex: contract, et1h')


    if args.build_b:
        build_basis(graph_complex)
    if args.build_op:
        build_operator(graph_complex)
    if args.rank:
        rank(graph_complex)
    if args.cohomology:
        cohomology(graph_complex)
    if args.square_zero:
        square_zero_test(graph_complex)
    if args.anti_commute:
        test_anti_commutativity(graph_complex)
    if args.commute:
        test_commutativity(graph_complex)
