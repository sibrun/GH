"""Compute the cohomology dimensions of graph complexes and bicomplexes.

Build graph complexes and bicomplexes.
Test whether an operator squares to zero, i.e. is a differential.
Test whether two operators anti-commute or commute, i.e. build a differential for a bicomplex.
Plot the cohomology dimensions of a graph complex.

This module parses command line arguments and contains the main function.
Building the basis and operator matrix as well as computing the matrix rank can be done in parallel for
different parameters. Use the option (-n_jobs) to specify the number of parallel jobs.
However, building a single basis or matrix file or computing a rank for a specific matrix can not be done in parallel.

There are options to ignore existing files (-ignore_ex), to display information (-info), plot information to a html file
(-plot_info), to show a progress bar (-pbar), for logging (-log warning), and for profiling (-profile).

Display information and show a progress bar are only available if not several processes work in parallel
(default: -n_jobs=1).

There are different options to compute the rank of operator matrices.
Use sage to determine the exact rank over the integers (-sage integer) or over a finite field (-sage mod).
Use linbox to determine the exact rank over the rationals (-linbox rational) or over a finite field (-linbox mod).
Use rheinfall to determine the exact rank over Z (-rheinfall mpz), over the rationals (-rheinfall mpq) or modulo 64
bit integers (-rheinfall int64).

There

Prerequisites:
    sagemath: This software is based on the sage math library.
        Download sagemath from:

            http://www.sagemath.org

        Activate a sage shell with

            $ sage -sh

            or use the command $ sage --python ... instead of $ python ... .

    tqdm: Install the tqdm module in the sage python version:

            (sage-sh)$ pip -install tqdm

    pandas: Install the pandas module in the sage python version:

            (sage-sh)$ pip -install pandas

    nauty: Used to generate graphs up to isomorphisms.
        Download nauty from:

            http://pallini.di.uniroma1.it/

    linbox: Used to determine the rank of matrices.
        Download linbox from:

            https://github.com/linbox-team

    rheinfall: Used to determine the rank of matrices.
        Download rheinfall from:

            https://github.com/riccardomurri/rheinfall

    For rank computations the corresponding example code of linbox or rheinfall is used.
    Replace the rank files by the ones provided in the folder GH/linbox_rheinfall_rank and compile them.
    Finally, move the executables to the folder GH/rank_exe.


Examples:
    In order to run the Graph Homology code, change to the directory GH, activate sage, and run:

    Ordinary graph complex:
        Build the basis:

            $ python ./source/GraphHomology.py ordinary -op1 contract -v 3,12 -l 0,9 -odd_e -n_jobs 4 -build_b

        Build the operator matrices for the differentials 'contract edges' and 'delete edges':

            $ python ./source/GraphHomology.py ordinary -op1 contract -v 3,12 -l 0,9 -odd_e -n_jobs 4 -build_op
            $ python ./source/GraphHomology.py ordinary -op1 delete -v 3,12 -l 0,9 -odd_e -n_jobs 4 -build_op

        Compute the ranks of the operator matrices:

            $ python ./source/GraphHomology.py ordinary -op1 contract -v 3,12 -l 0,9 -odd_e -n_jobs 4 -sage mod -rank
            $ python ./source/GraphHomology.py ordinary -op1 delete -v 3,12 -l 0,9 -odd_e -n_jobs 4 -sage mod -rank

        Test whether the operators square to zero, i.e. build a differential, test whether the operators anti-commute, and
            plot the cohomology dimensions of the respective graph complexes:

            $ python ./source/GraphHomology.py ordinary -op1 contract -op2 delete -v 0,12 -l 0,9 -odd_e -square_zero -anti_commute -cohomology

    Ordinary graph bicomplex 'contract_delete' with the differentials 'contract edges' and 'delete edges':

        $ python ./source/GraphHomology.py ordinary -bicomplex contract_delete -d 0,18 -odd_e -n_jobs 4 -sage mod -build_b -build_op -rank -square_zero -cohomology -csv

    Hairy graph complex:
        Build the basis:

            $ python ./source/GraphHomology.py hairy -op1 contract -v 3,11 -l 0,10 -hairs 0,9 -odd_e -odd_h -n_jobs 4 -build_b

        Build the operator matrices for the differentials 'contract edges' and 'edge to one hair':

            $ python ./source/GraphHomology.py hairy -op1 contract -v 3,11 -l 0,10 -hairs 0,9 -odd_e -odd_h -n_jobs 4 -build_op
            $ python ./source/GraphHomology.py hairy -op1 et1h -v 3,11 -l 0,10 -hairs 0,9 -odd_e -odd_h -n_jobs 4 -build_op

        Compute the ranks of the operator matrices:

            $ python ./source/GraphHomology.py hairy -op1 contract -v 3,11 -l 0,10 -hairs 0,9 -odd_e -odd_h -n_jobs 4 -sage mod -rank
            $ python ./source/GraphHomology.py hairy -op1 et1h -v 3,11 -l 0,10 -hairs 0,9 -odd_e -odd_h -n_jobs 4 -sage mod -rank

        Test whether the operators square to zero, i.e. build a differential, test whether the operators anti-commute, and
            plot the cohomology dimensions of the respective graph complexes:

            $ python ./source/GraphHomology.py hairy -op1 contract -op2 et1h -v 3,11 -l 0,10 -hairs 0,9 -odd_e -odd_h -square_zero -anti_commute -cohomology

    Hairy graph bicomplex 'contract_et1h' with the differentials 'contract edges' and 'edge to one hair':

        $ python ./source/GraphHomology.py hairy -bicomplex contract_et1h -d 0,14 -h_min ,-12,-1 -odd_e -odd_h -n_jobs 4 -sage mod -build_b -build_op -rank -square_zero -cohomology

    Bi colored hairy graph complex:
        Build the basis:

            $ python ./source/GraphHomology.py bi_c_hairy -op1 contract -v 3,8 -l 0,5 -hairs_a 0,6 -hairs_b 0,6 -odd_e -even_h_a -even_h_b -n_jobs 4 -build_b

        Build the operator matrices for the differentials 'contract edges' and 'split edges':

            $ python ./source/GraphHomology.py bi_c_hairy -op1 contract -v 3,8 -l 0,5 -hairs_a 0,6 -hairs_b 0,6 -odd_e -even_h_a -even_h_b -n_jobs 4 -build_op
            $ python ./source/GraphHomology.py bi_c_hairy -op1 split -v 3,8 -l 0,5 -hairs_a 0,6 -hairs_b 0,6 -odd_e -even_h_a -even_h_b -n_jobs 4 -build_op

        Compute the ranks of the operator matrices:

            $ python ./source/GraphHomology.py bi_c_hairy -op1 contract -v 3,8 -l 0,5 -hairs_a 0,6 -hairs_b 0,6 -odd_e -even_h_a -even_h_b -n_jobs 4 -sage mod -rank
            $ python ./source/GraphHomology.py bi_c_hairy -op1 split -v 3,8 -l 0,5 -hairs_a 0,6 -hairs_b 0,6 -odd_e -even_h_a -even_h_b -n_jobs 4 -sage mod -rank

        Test whether the operators square to zero, i.e. build a differential, test whether the operators anti-commute, and
            plot the cohomology dimensions of the respective graph complexes:

            $ python ./source/GraphHomology.py bi_c_hairy -op1 contract -op2 split -v 3,8 -l 0,5 -hairs_a 0,6 -hairs_b 0,6 -odd_e -even_h_a -even_h_b -square_zero -anti_commute -cohomology

    Bi colored hairy graph bicomplex 'contract_split' with the differentials 'contract edges' and 'split edges':

        $ python ./source/GraphHomology.py bi_c_hairy -bicomplex contract_split -d 0,11 -h_a_min ,-9,1 -h_b_min ,-9,1 -odd_e -even_h_a -even_h_b -n_jobs 4 -sage mod -build_b -build_op -rank -square_zero -cohomology

    Forested Graph complex:
        Build basis and operators:

        $ sage --python source/GraphHomology.py forested -v 0,7 -l 4,5 -marked 0,6 -hairs 0,1 -even_e -op1 contract -op2 unmark -build_b -build_op
"""

import argparse
import Log
import Profiling
import OrdinaryGraphComplex
import OrdinaryGraphBiComplex
import HairyGraphComplex
import HairyGraphBiComplex
import BiColoredHairyGraphComplex
import BiColoredHairyGraphBiComplex
import Parameters
import LinboxInterface
import RheinfallInterface
import CHairyGraphComplex
import ForestedGraphComplex
import WRHairyGraphComplex
import OrdinaryMerkulovComplex


logger = Log.logger.getChild('main')


def positive_int(value):
    value = int(value)
    if value <= 0:
        raise argparse.ArgumentTypeError('positive integer expected')
    return value


def int_value(arg):
    temp = arg.split(',')
    if len(temp) == 2:
        (prefix, value) = temp
    else:
        raise argparse.ArgumentTypeError('integer ,[-]value expected')
    return int(value)


def non_negative_range_type(arg):
    temp = arg.split(',')
    if len(temp) != 2:
        raise argparse.ArgumentTypeError(
            'range min,max with 0 < min < max expected')
    (min, max) = map(int, temp)
    if min >= max:
        raise argparse.ArgumentTypeError(
            'range min,max with 0 < min < max expected')
    return range(min, max)


def range_type(arg):
    temp = arg.split(',')
    if len(temp) == 3:
        (prefix, min, max) = temp
    elif len(temp) == 2:
        (min, max) = temp
    else:
        raise argparse.ArgumentTypeError(
            'range [,-]min,max with min < max expected')
    (min, max) = map(int, (min, max))
    if min >= max:
        raise argparse.ArgumentTypeError(
            'range [,-]min,max with min < max expected')
    return range(min, max)


def linbox_options(arg):
    options = [str(item) for item in arg.split(',')]
    if not set(options) <= LinboxInterface.linbox_options:
        raise argparse.ArgumentTypeError(
            'linbox options: ' + str(LinboxInterface.linbox_options))
    return options


def rheinfall_options(arg):
    options = [str(item) for item in arg.split(',')]
    if not set(options) <= RheinfallInterface.rheinfall_options:
        raise argparse.ArgumentTypeError(
            'rheinfall options: ' + str(RheinfallInterface.rheinfall_options))
    return options


def sage_rank_options(arg):
    options = [str(item) for item in arg.split(',')]
    if not set(options) <= Parameters.sage_rank_options:
        raise argparse.ArgumentTypeError(
            'sage rank options: ' + str(Parameters.sage_rank_options))
    return options


graph_types = ['ordinary', 'hairy', 'bi_c_hairy',
               'forested', 'wrhairy', 'chairy', 'ordinaryme']
operators = ['contract', 'delete', 'et1h', 'split',
             'unmark', 'contract_iso', 'unmark_iso']
bicomplexes = ['contract_et1h', 'contract_delete',
               'contract_split', 'contract_unmark', 'contract_unmark_iso', 'contract_unmark_top', ]

parser = argparse.ArgumentParser(
    description='Compute the homology of a graph complex')

parser.add_argument('graph_type', type=str,
                    choices=graph_types, help='type of the graphs')
parser.add_argument('-bicomplex', type=str,
                    choices=bicomplexes, help='bicomplex')
parser.add_argument('-op1', type=str, choices=operators, help='operator 1')
parser.add_argument('-op2', type=str, choices=operators,
                    default=None, help='operator 2')
parser.add_argument('-even_e', action='store_true', help='even edges')
parser.add_argument('-odd_e', action='store_true', help='odd edges')
parser.add_argument('-even_h', action='store_true', help='even hairs')
parser.add_argument('-odd_h', action='store_true', help='odd hairs')
parser.add_argument('-even_h_a', action='store_true', help='even hairs_a')
parser.add_argument('-odd_h_a', action='store_true', help='odd hairs_a')
parser.add_argument('-even_h_b', action='store_true', help='even hairs_b')
parser.add_argument('-odd_h_b', action='store_true', help='odd hairs_b')
parser.add_argument('-v', type=non_negative_range_type,
                    help='range min,max for number of vertices')
parser.add_argument('-l', type=non_negative_range_type,
                    help='range min,max for number of loops')
parser.add_argument('-marked', type=non_negative_range_type,
                    help='range min,max for number of marked edges')
parser.add_argument('-omega', type=non_negative_range_type,
                    help='range min,max for number of omega vertices')
parser.add_argument('-shift', type=int, default=1,
                    help='maximal shift = loops - vertices')
parser.add_argument('-hairs', type=non_negative_range_type,
                    help='range min,max for number of hairs')
parser.add_argument('-hairs_a', type=non_negative_range_type,
                    help='range min,max for number of hairs_a')
parser.add_argument('-hairs_b', type=non_negative_range_type,
                    help='range min,max for number of hairs_b')
parser.add_argument('-d', type=non_negative_range_type,
                    help='range min,max for degree of degree slices in bicomplex')
parser.add_argument('-h_min', type=range_type,
                    help='range min,max for minimal number of hairs in a degree slice of a bicomplex')
parser.add_argument('-h_a_min', type=range_type,
                    help='range min,max for minimal number of hairs_a in a degree slice of a bicomplex')
parser.add_argument('-h_b_min', type=range_type,
                    help='range min,max for minimal number of hairs_b in a degree slice of a bicomplex')
parser.add_argument('-ignore_ex', action='store_true',
                    help='ignore existing files')
parser.add_argument('-n_jobs', type=int, default=1,
                    help='number of parallel processes')
parser.add_argument('-pbar', action='store_true', help='show progressbar')
parser.add_argument('-profile', action='store_true', help='profiling')
parser.add_argument(
    '-log', type=str, choices=Log.log_levels_dict.keys(), help='logging level')
parser.add_argument('-info', action='store_true',
                    help='display info during calculations in browser')
parser.add_argument('-sage', type=sage_rank_options,
                    help='compute matrix ranks using the sage library, options: ' + str(Parameters.sage_rank_options))
parser.add_argument('-linbox', type=linbox_options,
                    help='compute matrix ranks using the linbox library, options: ' + str(LinboxInterface.linbox_options))
parser.add_argument('-rheinfall', type=str, choices=RheinfallInterface.rheinfall_options,
                    help="compute matrix ranks using the rheinfall library, options: " + str(RheinfallInterface.rheinfall_options))
parser.add_argument('-build', action='store_true',
                    help='build vector space basis and operator matrix')
parser.add_argument('-build_b', action='store_true',
                    help='build vector space basis')
parser.add_argument('-build_op', action='store_true',
                    help='build operator matrix')
parser.add_argument('-rank', action='store_true', help='compute matrix ranks')
parser.add_argument('-cohomology', action='store_true',
                    help='compute cohomology dimensions')
parser.add_argument('-csv', action='store_true',
                    help='export cohomology dimension to csv file')
parser.add_argument('-html', action='store_true',
                    help='export cohomology dimension to html file')
parser.add_argument('-square_zero', action='store_true',
                    help='square zero test')
parser.add_argument('-anti_commute', action='store_true',
                    help='test anti-commutativity of differentials')
parser.add_argument('-commute', action='store_true',
                    help='test commutativity of differentials')
parser.add_argument('-plot_info', action='store_true',
                    help='plot information about vector spaces and operator matrices')

args = parser.parse_args()


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def build_basis(graph_complex):
    logger.warning("\n----- Build Basis -----\n")
    graph_complex.build_basis(ignore_existing_files=args.ignore_ex, n_jobs=args.n_jobs, progress_bar=args.pbar,
                              info_tracker=args.info)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def build_operator(graph_complex):
    logger.warning("\n----- Build Matrix -----\n")
    graph_complex.build_matrix(ignore_existing_files=args.ignore_ex, n_jobs=args.n_jobs, progress_bar=args.pbar,
                               info_tracker=args.info)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def square_zero_test(graph_complex):
    logger.warning("\n----- Square Zero Test -----\n")
    graph_complex.square_zero_test()


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def test_anti_commutativity(graph_complex):
    logger.warning("\n----- Anti-commutativity Test -----\n")
    graph_complex.test_pairwise_anti_commutativity(commute=False)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def test_commutativity(graph_complex):
    logger.warning("\n----- Commutativity Test -----\n")
    graph_complex.test_pairwise_anti_commutativity(commute=True)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def rank(graph_complex):
    logger.warning("\n----- Compute Ranks -----\n")
    graph_complex.compute_rank(sage=args.sage, linbox=args.linbox, rheinfall=args.rheinfall,
                               ignore_existing_files=args.ignore_ex, n_jobs=args.n_jobs, info_tracker=args.info)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def cohomology(graph_complex):
    logger.warning("\n----- Compute Cohomology -----\n")
    graph_complex.plot_cohomology_dim(to_csv=args.csv, to_html=args.html)


@Profiling.cond_decorator(args.profile, Profiling.profile(Parameters.log_dir))
def plot_info(graph_complex):
    logger.warning("\n----- Plot Info -----\n")
    graph_complex.plot_info()


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

    logger.warning("\n###########################\n" +
                   "----- Graph Homology -----")

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
            if args.bicomplex == 'contract_delete':
                if args.d is None:
                    raise MissingArgumentError(
                        'specify -d: range for degree of degree slices in bicomplex')

                graph_complex = OrdinaryGraphBiComplex.OrdinaryContractDeleteBiGC(
                    args.d, even_edges)
            else:
                raise ValueError('Ordinary graphs bicomplex: contract_delete')

        elif len(operators) > 0 and set(operators) <= {'contract', 'delete'}:
            if args.v is None:
                raise MissingArgumentError(
                    'specify -v: range for number of vertices')
            if args.l is None:
                raise MissingArgumentError(
                    'specify -l: range for number of loops')

            graph_complex = OrdinaryGraphComplex.OrdinaryGC(args.v, args.l, even_edges, operators,
                                                            shift_loops_minus_vertices=args.shift)
        else:
            raise ValueError(
                'Differentials for ordinary graph complex: contract, delete')

    elif args.graph_type == 'ordinaryme':
        if len(operators) > 0 and set(operators) <= {'contract'}:
            if args.v is None:
                raise MissingArgumentError(
                    'specify -v: range for number of vertices')
            if args.l is None:
                raise MissingArgumentError(
                    'specify -l: range for number of loops')

            graph_complex = OrdinaryMerkulovComplex.OrdinaryMerkulovGC(args.v, args.l, even_edges, operators)
        else:
            raise ValueError(
                'Differentials for ordinary Merkulov graph complex: contract')

    elif args.graph_type == 'hairy':
        if args.even_h:
            even_hairs = True
        elif args.odd_h:
            even_hairs = False
        else:
            raise MissingArgumentError('specify -even_h or -odd_h')

        if args.bicomplex is not None:
            if args.bicomplex == 'contract_et1h':
                if args.d is None:
                    raise MissingArgumentError(
                        'specify -d: range for degree of degree slices in bicomplex')
                if args.h_min is None:
                    raise MissingArgumentError(
                        'specify -h_min: range for minimal number of hairs of degree slices in bicomplex')

                graph_complex = HairyGraphBiComplex.HairyCeEt1hBiGC(
                    args.d, args.h_min, even_edges, even_hairs)
            else:
                raise ValueError('Hairy graphs bicomplexes: ce_et1h')

        elif len(operators) > 0 and set(operators) <= {'contract', 'et1h'}:
            if args.v is None:
                raise MissingArgumentError(
                    'specify -v: range for number of vertices')
            if args.l is None:
                raise MissingArgumentError(
                    'specify -l: range for number of loops')
            if args.hairs is None:
                raise MissingArgumentError(
                    'specify -hairs: range for number of hairs')

            graph_complex = HairyGraphComplex.HairyGC(
                args.v, args.l, args.hairs, even_edges, even_hairs, operators)
        else:
            raise ValueError(
                'Differentials for hairy graph complex: contract, et1h')

    elif args.graph_type == 'bi_c_hairy':
        if args.even_h_a:
            even_hairs_a = True
        elif args.odd_h_a:
            even_hairs_a = False
        else:
            raise MissingArgumentError('specify -even_h_a or -odd_h_a')

        if args.even_h_b:
            even_hairs_b = True
        elif args.odd_h_b:
            even_hairs_b = False
        else:
            raise MissingArgumentError('specify -even_h_b or -odd_h_b')

        if args.bicomplex is not None:
            if args.bicomplex == 'contract_split':
                if args.d is None:
                    raise MissingArgumentError(
                        'specify -d: range for degree of degree slices in bicomplex')
                if args.h_a_min is None:
                    raise MissingArgumentError(
                        'specify -h_a_min: range for minimal number of hairs_a of degree slices in bicomplex')
                if args.h_b_min is None:
                    raise MissingArgumentError(
                        'specify -h_b_min: range for minimal number of hairs_b of degree slices in bicomplex')

                graph_complex = BiColoredHairyGraphBiComplex.BiColoredHairyContractSplitBiGC(args.d, args.h_a_min,
                                                                                             args.h_b_min, even_edges,
                                                                                             even_hairs_a, even_hairs_b)
            else:
                raise ValueError('Hairy graphs bicomplexes: ce_et1h')

        elif len(operators) > 0 and set(operators) <= {'contract', 'split'}:
            if args.v is None:
                raise MissingArgumentError(
                    'specify -v: range for number of vertices')
            if args.l is None:
                raise MissingArgumentError(
                    'specify -l: range for number of loops')
            if args.hairs_a is None:
                raise MissingArgumentError(
                    'specify -hairs_a: range for number of hairs_a')
            if args.hairs_b is None:
                raise MissingArgumentError(
                    'specify -hairs_b: range for number of hairs_b')

            graph_complex = BiColoredHairyGraphComplex.BiColoredHairyGC(args.v, args.l, args.hairs_a, args.hairs_b,
                                                                        even_edges, even_hairs_a, even_hairs_b, operators)
        else:
            raise ValueError(
                'Differentials for hairy graph complex: contract, split')

    elif args.graph_type == 'chairy':
        if len(operators) > 0 and set(operators) <= {'contract', 'contract_iso'}:
            if args.v is None:
                raise MissingArgumentError(
                    'specify -v: range for number of vertices')
            if args.l is None:
                raise MissingArgumentError(
                    'specify -l: range for number of loops')
            if args.hairs is None:
                raise MissingArgumentError(
                    'specify -hairs: range for number of hairs')

            graph_complex = CHairyGraphComplex.CHairyGC(
                args.v, args.l, args.hairs, even_edges, operators)

        else:
            raise ValueError(
                'Differentials for chairy graph complex: contract')
    elif args.graph_type == 'forested':
        if args.bicomplex is not None:
            if args.bicomplex == 'contract_unmark' or args.bicomplex == 'contract_unmark_iso':
                if args.l is None:
                    raise MissingArgumentError(
                        'specify -l: range for number of loops')
                if args.marked is None:
                    raise MissingArgumentError(
                        'specify -marked: range for number of marked edges')
                if args.hairs is None:
                    raise MissingArgumentError(
                        'specify -hairs: range for number of hairs')

                graph_complex = ForestedGraphComplex.ForestedContractUnmarkBiGC(
                    args.l, args.marked, args.hairs, even_edges, args.bicomplex == 'contract_unmark_iso')
            elif args.bicomplex == 'contract_unmark_top':
                if args.l is None:
                    raise MissingArgumentError(
                        'specify -l: range for number of loops')
                if args.marked is None:
                    raise MissingArgumentError(
                        'specify -marked: range for number of marked edges')
                if args.hairs is None:
                    raise MissingArgumentError(
                        'specify -hairs: range for number of hairs')
                graph_complex = ForestedGraphComplex.ContractUnmarkTopD(
                    args.l, args.marked, args.hairs, even_edges)
            else:
                raise ValueError(
                    'forested graphs bicomplex: contract_unmark or contract_unmark_iso')
        elif len(operators) > 0 and set(operators) <= {'contract', 'unmark'}:
            if args.v is None:
                raise MissingArgumentError(
                    'specify -v: range for number of vertices')
            if args.l is None:
                raise MissingArgumentError(
                    'specify -l: range for number of loops')
            if args.marked is None:
                raise MissingArgumentError(
                    'specify -marked: range for number of marked edges')
            if args.hairs is None:
                raise MissingArgumentError(
                    'specify -hairs: range for number of hairs')

            graph_complex = ForestedGraphComplex.ForestedGC(
                args.v, args.l, args.marked, args.hairs, even_edges, operators)

        else:
            raise ValueError(
                'Differentials for forested graph complex: contract, unmark')
    elif args.graph_type == 'wrhairy':
        if len(operators) > 0 and set(operators) <= {'contract', 'contract_iso'}:
            if args.v is None:
                raise MissingArgumentError(
                    'specify -v: range for number of vertices')
            if args.l is None:
                raise MissingArgumentError(
                    'specify -l: range for number of loops')
            if args.hairs is None:
                raise MissingArgumentError(
                    'specify -hairs: range for number of hairs')
            if args.omega is None:
                raise MissingArgumentError(
                    'specify -omega: range for number of omega vertices')

            graph_complex = WRHairyGraphComplex.WRHairyGC(
                args.v, args.l, args.hairs, args.omega, operators)

        else:
            raise ValueError(
                'Differentials for wrhairy graph complex: contract')

    if args.plot_info:
        plot_info(graph_complex)
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
