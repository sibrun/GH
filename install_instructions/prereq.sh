#! /bin/sh
#
## Customization section

# comment any of the following lines if you already have
# the required version (or a newer one) installed in the
# standard PATH/LD_LIBRARY_PATH
#
#CYTHON=0.11.1
#PYTHON=2.5.4
#SWIG=1.3.35
FFLAS_FFPACK=2.2.0 # https://github.com/linbox-team/fflas-ffpack
GMP=6.1.0
GIVARO=4.0.1 # https://github.com/linbox-team/givaro/archive/v4.0.1.tar.gz
LINBOX=1.4.0
#ATLAS=3.8.3
BOOST=1.60.0 # https://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.tar.bz2/download#
#TCMALLOC=2.5 # https://github.com/gperftools/gperftools
#TBB=40_278oss # http://threadingbuildingblocks.org/uploads/78/179/4.0%20update%202/tbb40_278oss_lin.tgz

# set this to install BeBOP's Sparse Matrix Converter library and command-line utility
#SMC=yes # See: http://bebop.cs.berkeley.edu/smc/

# if your cluster is homogeneous (i.e., nodes all have the same
# architecture and memory size), undefining this might give some speedup
#BOOST_MPI_HOMOGENEOUS=yes


## No customization should be necessary further down here

PROG="$(basename $0)"

usage () {
cat <<EOF
Usage: $PROG [DIR [CFLAGS ...]]
Download and install all required software for 'rheinfall' into
directory DIR.  If omitted, DIR defaults to '`pwd`/sw'.
Any additional arguments are placed directly on the 'configure'
command line of every package, so it can be used to set CFLAGS,
CXXFLAGS etc.
EOF
}


## helper functions
die () {
  rc="$1"
  shift
  (echo -n "$PROG: ERROR: ";
      if [ $# -gt 0 ]; then echo "$@"; else cat; fi) 1>&2
  exit $rc
}

have_command () {
  type "$1" >/dev/null 2>/dev/null
}

require_command () {
  if ! have_command "$1"; then
    die 1 "Could not find required command '$1' in system PATH. Aborting."
  fi
}

is_absolute_path () {
  expr "$1" : '/' >/dev/null 2>/dev/null
}

_ () {
    echo
    echo ==== "$@" ...;
}


## parse command-line

case "$1" in
    -h|--help)
        usage
        exit 0
        ;;
    *)
esac


## main

require_command bzip2
require_command gzip
require_command make
require_command tar
require_command wget

set -e

# find Rheinfall source directory
if [ -e util/${PROG} ]; then
    etc_dir=$(pwd)/util
elif [ -e ${PROG} ]; then
    etc_dir=$(pwd)
fi

# target directory
if [ -n "$1" ]; then
    sw="$1"
    shift
else
    sw='sw'
fi
if ! is_absolute_path "$sw"; then
    sw="$(pwd)/${sw}"
fi
mkdir -p "${sw}"
mkdir -p "${sw}/src"
cd "${sw}"

# paths
PATH=${sw}/bin:$PATH
LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH
PYTHONPATH=${sw}/lib/python

# misc
ncpus="$(grep -c '^processor' /proc/cpuinfo)"
concurrent_make="make -j $ncpus"

# the `cpu MHz` line in /proc/cpuinfo states the *current* CPU clock,
# whereas "bogomips" is approximately twice the maximal clock on all Intel/AMD CPUs,
# see: http://tldp.org/HOWTO/BogoMips/bogo-faq.html#AEN192
bogomips=`fgrep bogomips /proc/cpuinfo | head -1 | cut -d: -f2 | cut -d. -f1`
mhz=`expr $bogomips / 2`

case `uname -m` in
        i?86) bits=32 ;;
        x86_64) bits=64 ;;
        *) die 1 "Unknown architecture `uname -m`: is it 32-bit or 64-bit?" ;;
esac

if egrep '^flags: .*\<avx\>' /proc/cpuinfo; then
  maybe_enable_avx='--enable-avx'
fi

if egrep '^flags: .*\<sse4\>' /proc/cpuinfo; then
  maybe_enable_sse='--enable-sse'
fi

# Python
if [ -n "${PYTHON}" ]; then
    _ Installing Python ${PYTHON}
    cd ${sw}/src/
    wget -N http://www.python.org/ftp/python/${PYTHON}/Python-${PYTHON}.tar.bz2
    set -x
    tar -xjf Python-${PYTHON}.tar.bz2
    cd Python-${PYTHON}
    ./configure --prefix=${sw} "$@"
    $concurrent_make
    make install
    set +x
fi # PYTHON


# SWIG
if [ -n "${SWIG}" ]; then
    _ Downloading SWIG ${SWIG}
    cd ${sw}/src/
    wget -N http://mirror.switch.ch/ftp/ubuntu/pool/main/s/swig1.3/swig1.3_${SWIG}.orig.tar.gz
    set -x
    tar -xzf swig1.3_${SWIG}.orig.tar.gz
    cd swig-${SWIG}
    ./configure --prefix=${sw} \
        --with-python=${sw}/bin \
        --without-allegrocl \
        --without-chicken \
        --without-clisp \
        --without-csharp \
        --without-gcj \
        --without-guile \
        --without-java \
        --without-lua \
        --without-mzscheme \
        --without-ocaml \
        --without-octave \
        --without-perl5 \
        --without-php4 \
        --without-pike \
        --without-r \
        --without-ruby \
        --without-rxspencer \
        --without-tcl \
        "$@";
    $concurrent_make
    make install
    set +x
fi # SWIG


# Cython
if [ -n "$CYTHON" ]; then
    _ Installing Cython-${CYTHON}
    cd ${sw}/src/
    wget -N http://www.cython.org/Cython-${CYTHON}.tar.gz
    set -x
    tar -xzf Cython-${CYTHON}.tar.gz
    cd Cython-${CYTHON}
    python setup.py build
    python setup.py install --home=${sw}

    PYTHONPATH=$PYTHONPATH:`pwd`; export PYTHONPATH
    PATH=$PATH:`pwd`/bin; export PATH
    set +x
fi # CYTHON


# GMP
if [ -n "$GMP" ]; then
    _ Installing GMP ${GMP}
    cd ${sw}/src
    #wget -N http://mirror.switch.ch/ftp/mirror/gnu/gmp/gmp-${GMP}.tar.bz2
    wget -N https://gmplib.org/download/gmp/gmp-${GMP}.tar.bz2
    set -x
    tar -xjf gmp-${GMP}.tar.bz2
    cd gmp-${GMP}
    ./configure --prefix=${sw} --enable-cxx "$@"
    $concurrent_make
    make install
    set +x
fi # GMP


# GIVARO (cfr. http://groups.google.com/group/linbox-use/browse_thread/thread/82673844f6921271)
if [ -n "$GIVARO" ]; then
    _ Installing Givaro ${GIVARO}
    cd ${sw}/src
    wget -N https://github.com/linbox-team/givaro/archive/v${GIVARO}.tar.gz -O givaro-${GIVARO}.tar.gz
    #old: wget -N http://www-lmc.imag.fr/CASYS/LOGICIELS/givaro/Downloads/givaro-${GIVARO}.tar.gz
    set -x
    tar -xzf givaro-${GIVARO}.tar.gz
    cd givaro-${GIVARO%%rc[0-9]}
    # work around bug in ./configure: the test for GMP cannot find it
    # unless it's in the LD_LIBRARY_PATH
    export LD_LIBRARY_PATH=${sw}/lib:$LD_LIBRARY_PATH
    ./autogen.sh
    ./configure  --prefix=${sw} --enable-shared ${GMP:+"--with-gmp=${sw}"} "$@"
    $concurrent_make
    make install
    set +x
fi # GIVARO


# FFLAS-FFPACK
if [ -n "$FFLAS_FFPACK" ]; then
    _ Installing FFLAS-FFPACK ${FFLAS_FFPACK}
    cd ${sw}/src
    if ! [ -d fflas-ffpack-${FFLAS_FFPACK} ]; then
      git clone https://github.com/linbox-team/fflas-ffpack.git fflas-ffpack-${FFLAS_FFPACK}
    fi
    set -x
    cd fflas-ffpack-${FFLAS_FFPACK}
    git checkout v${FFLAS_FFPACK}
    ./autogen.sh --prefix=${sw} \
	${GIVARO:+"--with-givaro=${sw}"} \
	--with-blas-libs="-lopenblas" \
        --enable-openmp \
	$maybe_enable_sse \
	$maybe_enable_avx \
        "$@";
    $concurrent_make
    make install
    set +x
fi # FFLAS-FFPACK


# ATLAS
if [ -n "$ATLAS" ]; then
    cd ${sw}/src
    wget "http://switch.dl.sourceforge.net/sourceforge/math-atlas/atlas${ATLAS}.tar.bz2" \
         -O atlas${ATLAS}.tar.bz2
    set -x
    tar -xjf atlas${ATLAS}.tar.bz2
    cd ATLAS
    mkdir -p BLDdir
    cd BLDdir
    ../configure -v 2 \
        -b ${bits} \
        -m ${mhz} \
        -D c -DPentiumCPS=${mhz} \
        -Si cputhrchk 0 \
        -Ss pmake "$concurrent_make" \
        --prefix=${sw}
    make build
    (cd lib; make cshared cptshared && cp -a *.so ${sw}/lib)
    #make check
    #make time
    make install
    set +x
fi # ATLAS


# LinBox
if [ -n "$LINBOX" ]; then
    _ Installing LinBox ${LINBOX}
    cd ${sw}/src
    #wget -N http://linalg.org/linbox-${LINBOX}.tar.gz
    #tar -xzf linbox-${LINBOX}.tar.gz
    if ! [ -d linbox-${LINBOX} ]; then
      git clone https://github.com/linbox-team/linbox.git linbox-${LINBOX}
    fi
    set -x
    cd linbox-${LINBOX}
    git checkout v${LINBOX}
    ./autogen.sh --prefix=${sw} \
        ${ATLAS:+"--with-blas=${sw}/lib"} \
        ${GMP:+"--with-gmp=${sw}"} \
        ${GIVARO:+"--with-givaro=${sw}"} \
	${FFLAS_FFPACK:+"--with-fflas-ffpack=${sw}"} \
        --enable-openmp \
        "$@";
    $concurrent_make
    make install
    set +x
fi # LINBOX


# Boost
if [ -n "$BOOST" ]; then
    _ Installing Boost ${BOOST}
    cd ${sw}/src
    boost_file=$(echo boost_$BOOST | tr . _)
    wget \
        "https://sourceforge.net/projects/boost/files/boost/${BOOST}/${boost_file}.tar.bz2/download#" \
        -O "${boost_file}.tar.bz2"
    set -x
    tar -xaf  "${boost_file}.tar.bz2"
    cd ${boost_file}
    #patch -p1 -i ${etc_dir}/boost_1_45_0.ssend.patch
    # build Boost.MPI for homogeneous clusters (same arch, so avoid pack/unpack)
    if [ "x$BOOST_MPI_HOMOGENEOUS" = "xyes" ]; then
        sed -e 's|^//#define BOOST_MPI_HOMOGENEOUS|#define BOOST_MPI_HOMOGENEOUS|' \
            -i boost/mpi/config.hpp
    fi
    ./bootstrap.sh --prefix=${sw} --with-libraries=mpi,serialization,test \
        variant=release threading=multi
    cat >> project-config.jam <<EOF
# Boost will not build Boost.MPI unless it is explicitly
# told to by the following line:
using mpi : mpicxx ;
EOF
    # build Boost with the new `bjam`
    PATH=$(pwd)/tools/jam/src/bin.$(uname -s | tr A-Z a-z)$(uname -m):$PATH
    export PATH
    ./bjam --prefix=${sw} threading=multi variant=release install
    set +x
fi # BOOST


# Google perftools
if [ -n "$TCMALLOC" ]; then
    _ Installing Google PerfTools $TCMALLOC ...
    cd ${sw}/src
    set -x
    git clone https://github.com/gperftools/gperftools.git
    #wget -N http://google-perftools.googlecode.com/files/google-perftools-${TCMALLOC}.tar.gz
    #tar -xzf "google-perftools-${TCMALLOC}.tar.gz"
    cd gperftools
    git checkout gperftools-${TCMALLOC}
    autoreconf -i
    ./configure --prefix=${sw} \
        --enable-frame-pointers --disable-debugalloc \
        "$@";
    $concurrent_make
    make install
    set +x
fi


# Intel TBB
if [ -n "$TBB" ]; then
    _ Installing Intel TBB $TBB ...
    cd ${sw}/src/
    set -x
    case "$TBB" in
        # TBB download URL changes with release ...
        40_278oss)
            wget -N "http://threadingbuildingblocks.org/uploads/78/179/4.0%20update%202/tbb40_278oss_lin.tgz"
            wget -N "http://threadingbuildingblocks.org/uploads/78/179/4.0%20update%202/tbb40_278oss_src.tgz"
            ;;
        30_20100915oss)
            wget -N "http://www.threadingbuildingblocks.org/uploads/77/161/3.0%20update%203/tbb30_20100915oss_lin.tgz"
            wget -N "http://www.threadingbuildingblocks.org/uploads/77/161/3.0%20update%203/tbb30_20100915oss_src.tgz"
            ;;
        *)
            die 1 "Unknown download URL for Intel TBB ${TBB}"
            ;;
    esac
    mkdir -p ${sw}/opt
    for tgz in "tbb${TBB}"*.tgz; do
        tar -xzf "$tgz" -C "${sw}/opt"
    done
    ln -sf "tbb${TBB}" "${sw}/opt/tbb"
    # compile
    tbb_root="${sw}/opt/tbb"
    cd "${tbb_root}"
    make clean default
    # there's no "install" make target, so we just copy files
    tbb_arch=$(make info | grep '^arch=' | cut -d= -f2)
    tbb_runtime=$(make info | grep '^runtime=' | cut -d= -f2)
    tbb_build_prefix=$(make info | grep '^tbb_build_prefix=' | cut -d= -f2)
    if [ -d "build/${tbb_build_prefix}_release" ]; then
        mkdir -p  "lib/${tbb_arch}/${tbb_runtime}/"
        cp "build/${tbb_build_prefix}_release"/lib* "lib/${tbb_arch}/${tbb_runtime}/"
    fi
    if [ -d "build/${tbb_build_prefix}_debug" ]; then
        mkdir -p "lib/${tbb_arch}/${tbb_runtime}/"
        cp "build/${tbb_build_prefix}_debug"/lib* "lib/${tbb_arch}/${tbb_runtime}/"
    fi
    set +x
fi


# BeBOP Sparse Matrix Converter
if [ -n "$SMC" ] && [ ! "$SMC" = no ]; then
    _ Installing BeBOP Sparse Matrix Converter ...
    cd ${sw}/src/
    set -x
    wget http://bebop.cs.berkeley.edu/smc/tarballs/bebop_make.tar.gz
    wget http://bebop.cs.berkeley.edu/smc/tarballs/bebop_util.tar.gz
    wget http://bebop.cs.berkeley.edu/smc/tarballs/sparse_matrix_converter.tar.gz
    for tgz in bebop_make.tar.gz bebop_util.tar.gz sparse_matrix_converter.tar.gz; do
        tar -xzf "$tgz"
    done
    # configure for running on Linux
    cd bebop_make
    rm Makefile.include; ln -s Makefile.include.linux Makefile.include
    #sed -i -e "s/^CC *= *gcc/CC = ${CC}/"
    #sed -i -e "s/^LINKER *= *gcc/LINKER = ${CC}/"
    cd ../bebop_util
    make
    cp -av libbebop_util.a  libbebop_util.so "${sw}/lib"
    mkdir -p "${sw}/include/bebop"
    cp -av include/bebop/util "${sw}/include/bebop/util"
    cd ../sparse_matrix_converter
    make
    cp -av sparse_matrix_converter "${sw}/bin"
    cp -av libsparse_matrix_converter.a libsparse_matrix_converter.so "${sw}/lib"
    cp -av include/bebop/smc "${sw}/include/bebop/smc"
fi

_ All done.