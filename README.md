feyntrop
==========================================

**feyntrop** is a computer program to evaluate Feynman integrals. The core `C++` integration code, written mainly by [Michael Borinsky](//michaelborinsky.com), is an update of the proof-of-concept implementation [tropical-feynman-quadrature](https://github.com/michibo/tropical-feynman-quadrature), published with the paper ['Tropical Monte Carlo quadrature for Feynman integrals'](//arxiv.org/abs/2008.12310). **feyntrop** can be used through a high-level [Python](//python.org) interface, written by Henrik Munch.

The whole **feyntrop** project  is managed  as [a github repository](//github.com/michibo/feyntrop/). Comments, bug reports, and pull requests are very welcome.

If **feyntrop** is helpful in your research, please cite [M. Borinsky, H. J. Munch, F. Tellander: 'Tropical Feynman integration in the Minkowski regime', *Computer Physics Communications*, 292 (2023), 108874](//doi.org/10.1016/j.cpc.2023.108874) [arXiv:2302.08955](//arxiv.org/abs/2302.08955)
as well as [M. Borinsky: 'Tropical Monte Carlo quadrature for Feynman integrals', Ann. Inst. Henri Poincaré Comb. Phys. Interact. 10 (2023), no. 4, pp. 635–685](//doi.org/10.4171/AIHPD/158) [arXiv:2008.12310](//arxiv.org/abs/2008.12310).

The implementation internally uses [Eigen](http://eigen.tuxfamily.org), OpenMP, the [xoshiro256+](http://prng.di.unimi.it/) random number generator, and the [JSON for Modern C++](//github.com/nlohmann/json) library.

Download
--------
 
To download **feyntrop** use the command

```
git clone https://github.com/michibo/feyntrop.git
cd feyntrop
```

Compilation
-----------

To (re)compile **feyntrop** call

```
make clean && make
```

If the compilation fails, check if a suitable `C++` compiler is installed (see below for macOS problems). 

### Problems specific to macOS

The default macOS `C++` compiler does not support OpenMP. So, on macOS, the compilation of **feyntrop** might also fail with an error mentioning the `-fopenmp` flag. If you use homebrew, installing an OpenMP-compatible `C++` compiler via the command
```
brew install libomp
```
might help. If that does not help, then the paths to the compiler binaries in the `Makefile` probably need to be adjusted for your local environment. See the `Makefile` for some hints on how to do this.

Tests
-----

To run the tests, make sure to have [Python](//python.org) installed.
To ensure that **feyntrop** was built correctly, run the file `/tests/test_suite.py`. That means running
```
cd tests
python test_suite.py
```
This Python script uses **feyntrop** to compute examples between 1-2 loops and 2-5 points and then compares them against pre-computed values.
Ratios between newly computed and pre-computed coefficients in the epsilon expansion will be printed, which should all be close to 1.

If you cannot or don't want to use `Python`, you can also directly test the `C++`-compiled code by running
```
./feyntrop < low_level_input.json
```
In the top-level directory. The output should be the value of the Feynman integral of the wheel graph with three spokes. This output roughly looks as follows (see the **Low-level interface** section below for details of the format):
```
{"IGtr":84.0,"integral":[[[7.215238614660525,0.00203586844683068],[0.0,0.0]],[[-57.629482716637696,0.018239280410844466],[0.0,0.0]],[[240.79344300578586,0.09697082078903732],[0.0,0.0]]],"seconds preprocessing":0.001007578,"seconds sampling":2.0064038810000002}
```

Tutorial
--------

To run the tutorial (which uses the high-level interface) in notebook form, you must install [jupyter-notebook](//jupyter.org/). Run
```
jupyter notebook tutorial_2L_3pt.ipynb
```
in the top directory of this repository to open the tutorial notebook.

Examples
--------

You might prefer a simple script file instead of a `jupyter notebook`. The file `simple_example_2L_3pt.py` contains the example from the tutorial notebook. The comments in this file completely explain the usage of `feyntrop` using the `Python` interface.

Moreover, the `examples/` folder in this repository contains a variety of more complicated examples for Feynman integral computations using `feyntrop` with the Python interface.  In this folder, you can, for instance, run the 5-loop 2-point example from the paper with the command 
```
python 5L_2pt.py
```

Low-level interface
-------------------

**feyntrop** can also be used with a low-level command-line interface without Python. This might be convenient in a high-performance computing environment. To use this interface, create a file similar to the `low_level_input.json` file in this repository. Here is the content of this file:
```
{
  "graph" : [ [ [0, 1], 1 ], [ [0, 2], 1 ], [ [0, 3], 1 ], [ [1, 2], 1 ], [ [2, 3], 1], [ [3, 1], 1 ] ],
  "dimension" : 4,
  "scalarproducts" : [ [ -3,  1,  1,  1 ],
                       [  1, -3,  1,  1 ],
                       [  1,  1, -3,  1 ],
                       [  1,  1,  1, -3 ] ],
  "masses_sqr" : [ 0, 1, 2, 3, 4, 5 ],
  "num_eps_terms" : 3,
  "lambda" : 0,
  "N" : 10000000,
  "seed" : 0
}
```
The field `"graph"` encodes the Feynman graph. It is a list of edges of the form 
```
[ [ [v0, w0], nu0 ], [ [v1, w1], nu1 ], [ [v2, w2], nu2 ],... ]
```
where `v0,w0` and so on are pairs of vertices corresponding to an edge and `nu0` is the corresponding edge weight.
The field `"scalarproducts"` is a matrix of scalar products. The `(v,w)`-th entry of the matrix is the scalar product of `p_u * p_v`, where `p_u` is the incoming momentum into vertex `u`. Hence, the matrix must be symmetric and have as many rows and columns as vertices. (Vertices without incoming momentum can be represented by setting the respective row and column equal to 0.) Due to momentum conservation, the rows and columns of the matrix must sum to 0.

The field `"masses_sqr"` is a list of masses containing one mass for each edge. (Of course, the masses might be 0.)
The field `"lambda"` is the deformation parameter, `"dimension"` is the spacetime dimension, `"num_eps_terms"` is the order in the epsilon expansion that should be computed, and `"N"` is the number of points that shall be sampled. The `"seed"` is the seed for the random number generators. For most practical purposes, the seed can be set to 0.

The content of the JSON file must be piped into the **feyntrop** executable file, which is created in the top directory of this repository by the `make` command. For instance, like this:
```
./feyntrop < low_level_input.json
```
Among some logging information (via stderr), this command produces the output (via stdout) in JSON format 
```
{"IGtr":84.0,"integral":[[[7.215238614660525,0.00203586844683068],[0.0,0.0]],[[-57.629482716637696,0.018239280410844466],[0.0,0.0]],[[240.79344300578586,0.09697082078903732],[0.0,0.0]]],"seconds preprocessing":0.001007578,"seconds sampling":2.0064038810000002}
```
where the field `"IGtr"` is the tropicalized Feynman integral (in this case equal to the *Hepp bound*), 
the field `"integral"` contains the epsilon expansion of the integral (without the usual gamma-function prefactor), `I = I0 + eps * I1  + eps^2 * I2 + ...`,
in the form
```
[ [ [ Re(I0), Delta(Re(I0)) ], [Im(I0), Delta(Im(I0))] ], [ [ Re(I1), Delta(Re(I1)) ], [Im(I1), Delta(Im(I1))] ], ... ]
```
where `Delta` is the respective error term (i.e. one expected standard deviation) and 
where `Re` and `Im` denote the real and imaginary part of the respective coefficient in the expansion.

The other fields store the sampling and the preprocessing time.

If you are not interested in the logging information, use the command
```
./feyntrop < low_level_input.json 2> /dev/null
```
instead.

Changing the number of CPUs used
--------------------------------

By default, **feyntrop** uses the maximal number of available CPUs in the sampling step. This behaviour can be changed using the environment variable `OMP_NUM_THREADS`.
For instance, the command
```
OMP_NUM_THREADS=2 ./feyntrop < low_level_input.json 2> /dev/null
```
performs the sampling step for integrating the input Feynman integral problem from `low_level_input.json` with only two threads.

A similar option can be used with the Python interface.  For instance,
```
OMP_NUM_THREADS=2 python 5L_2pt.py
```
runs the 5 loop 2-point integral example with two threads.
