# Graph Homolgy

Compute the cohomology dimensions of graph complexes.

## Getting Started

### Prerequisites
This graph homology library is based on the sage math library. Install the [sage math](http://www.sagemath.org) library and activate the sage shell:
```
$ sage -sh
```
Install the pandas and tqdm modules:
```
(sage-sh)$ pip -install tqdm
(sage-sh)$ pip -install pandas
```
Furthermore, in order to generate graphs up to isomorphism the graph library [nauty](http://pallini.di.uniroma1.it/) is needed.

Small matrices can be handled using the sage function to compute exact ranks.

In order to determine ranks of large matrices [linbox](https://github.com/linbox-team) and or [rheinfall](https://github.com/riccardomurri/rheinfall)
are required. 
For rank computations the corresponding example code of [linbox](https://github.com/linbox-team/linbox/blob/master/examples/rank.C) or 
[rheinfall](https://github.com/riccardomurri/rheinfall/blob/master/src.c%2B%2B/examples/rank.cpp) is used.
Replace the rank files by the ones provided in [this](https://github.com/sibrun/GH/tree/master/linbox_rheinfall_modified_rank) folder and compile them.
Finally move the executables to the folder [GH/source](https://github.com/sibrun/GH/tree/master/source).

## Running the tests

Change to the [GH](https://github.com/sibrun/GH) directuory to run the tests:
```
$ sage --python ./source/TestOrdinaryGraphComplex.py
$ sage --python ./source/TestOrdinaryGraphBiComplex.py
$ sage --python ./source/TestHairyGraphComplex.py
$ sage --python ./source/TestHairyGraphBiComplex.py
$ sage --python ./source/TestBiColoredHairyGraphComplex.py
$ sage --python ./source/TestBiColoredHairyGraphBiComplex.py
```

## Running the examples

The file [GraphHomology.py](https://github.com/sibrun/GH/blob/master/source/GraphHomology.py) provides commands to run the implemented
examples. 

Change to the directory GH and build the basis for the ordinary graph complex with odd edges:
```
$ sage --python ./source/GraphHomology.py ordinary -op1 contract -v 3,12 -l 0,9 -odd_e -n_jobs 4 -build_b
```
Afterwards build the operator matrices of the conctract edges and delete edges differentials: 
```
$ sage --python ./source/GraphHomology.py ordinary -op1 contract -v 3,12 -l 0,9 -odd_e -n_jobs 4 -build_op
$ sage --python ./source/GraphHomology.py ordinary -op1 delete -v 3,12 -l 0,9 -odd_e -n_jobs 4 -build_op
```
Then compute the ranks of the operator matrices:
```
$ sage --python ./source/GraphHomology.py ordinary -op1 contract -v 3,12 -l 0,9 -odd_e -n_jobs 4 -sage mod -rank
$ sage --python ./source/GraphHomology.py ordinary -op1 delete -v 3,12 -l 0,9 -odd_e -n_jobs 4 -sage mod -rank
```
Finally test whether the operators square to zero, i.e. build a differential, whether they anti-commute and plot the 
cohomology dimensions of the associated graph complexs for each differential:
```
$ sage --python ./source/GraphHomology.py ordinary -op1 contract -op2 delete -v 0,12 -l 0,9 -odd_e -square_zero -anti_commute -cohomology
```
The comhomology dimensions of the bicomplex with the two differentials contract edges and delete edges can be computed with:
```
$ sage --python ./source/GraphHomology.py ordinary -bicomplex contract_delete -d 0,18 -odd_e -n_jobs 4 -sage mod -build_b -build_op -rank -square_zero -cohomology -csv
```

More examples for hairy graph complexes can be found in the file [GraphHomology.py](https://github.com/sibrun/GH/blob/master/source/GraphHomology.py).

Feel free to extend this library with more examples of graph complexes.

## Documentation
In order to generate a documentation using [sphinx](http://www.sphinx-doc.org/en/master/#) change to the folder 
[GH/docs](https://github.com/sibrun/GH/tree/master/docs) activate the sage shell and generate a html documentation with:
```
(sage-sh)$ make html
```
The html documentation will be found under GH/docs/build/html/idex.html.

## Authors

* **Simon Brun** 
* **Thomas Willwacher**

## License


## Acknowledgments



