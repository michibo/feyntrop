[feyntrop](//michaelborinsky.com/feyntrop)
==========================================

**feyntrop** is a computer program to evaluate Feynman integrals. The core C++ integration code, written mainly by [Michael Borinsky](//michaelborinsky.com), is an update of the proof-of-concept implementation [tropical-feynman-quadrature](https://github.com/michibo/tropical-feynman-quadrature), which was published with the paper ['Tropical Monte Carlo quadrature for Feynman integrals'](//arxiv.org/abs/2008.12310). **feyntrop** can be used through a high-level [python](//python.org) interface, written by Henrik Munch.

The whole **feyntrop** program is available as [this github repository](//github.com/michibo/feyntrop/). Comments, bug-reports and pull-requests are very welcome.

If **feyntrop** is helpful in your research, please cite,
[M. Borinsky, H. J. Munch, F. Tellander: 'Tropical Feynman integration in the Minkowski regime' arXiv:2302.08955](//arxiv.org/abs/2302.08955) as well as [M. Borinsky: 'Tropical Monte Carlo quadrature for Feynman integrals' arXiv:2008.12310](//arxiv.org/abs/2008.12310).

The implementation internally uses [Eigen](http://eigen.tuxfamily.org), OpenMP, the [xoshiro256+](http://prng.di.unimi.it/) random number generator and the [JSON for Modern C++](//github.com/nlohmann/json) library.

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

If the compilation fails, check if you have a suitable C++ compiler installed (see below for macOS problems). 

### Problems specific to macOS

The default macOS C++ compiler does not support OpenMP. So on macOS, the compilation of **feyntrop** might also fail with an error mentioning the `-fopenmp` flag. If you use homebrew, installing an OpenMP compatible C++ compiler via the command
```
brew install libomp
```
might help. If that does not help, then probably the paths to the compiler binaries in the `Makefile` need to be adjusted for your local environment. See the `Makefile` for some hints on how to do this.

Tests
-----

To run the tests, make sure to have [python](//python.org) installed.
To ensure that **feyntrop** has been built correctly, please run the file `/tests/test_suite.py`. That means, run
```
cd tests
python test_suite.py
```
This python script uses **feyntrop** to compute examples between 1-2 loops and 2-5 points, and then compares against pre-computed values.
Ratios between newly computed and pre-computed coefficients in the epsilon expansion will be printed, which should all be close to 1.

Tutorial
--------

To run the tutorial (which uses the high-level interface) in notebook form, you have to have [jupyter-notebook](//jupyter.org/) installed. Run
```
jupyter notebook tutorial_2L_3pt.ipynb
```
in the top directory of this repository to open the tutorial notebook.

Examples
--------

Examples can be found in the `examples/` folder in this repository. 
In this folder, one can for instance run the 5 loop 2-point example from the paper with the command
```
python 5L_2pt.py
```

Low-level interface
-------------------

**feyntrop** can also be used without python. For instance, this might be convenient in a high-performance computing environment. To use this interface, create a file similar to the `low_level_input.json` file in this repository. Here is the content of this file:
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
  "N" : 10000000
}
```
The field `"graph"` encodes the Feynman graph. It is a list of edges of the form 
```
[ [ [v0, w0], nu0 ], [ [v1, w1], nu1 ], [ [v2, w2], nu2 ],... ]
```
where `v0,w0` and so on are pairs of vertices corresponding to an edge and `nu0` is the corresponding edge weight.
The field `"scalarproducts"` is a matrix of scalar products. The `(v,w)`-th entry of the matrix is the scalar product of `p_u * p_v` where `p_u` is the incoming momentum into vertex `u`. The matrix must hence be symmetric and have as many rows and columns as there are vertices. (Vertices without incoming momentum can be represented by setting the respective row and column equal to 0.) Due to momentum conservation, the rows and columns of the matrix must sum to 0.

The field `"masses_sqr"` is a list of masses, which contains one mass for each edge. (Of course, the masses might be 0.)
The field `"lambda"` is the deformation parameter, `"dimension"` is the spacetime dimension, `"num_eps_terms"` is the order in the epsilon expansion that should be computed and `"N"` is the number of points that shall be sampled.

The content of the JSON file must be piped into the **feyntrop** executable file, which is created in the top-directory of this repository by the `make` command. For instance, like this:
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

The other fields give store the sampling and the preprocessing time.

If you are not interested in the logging information, use, for instance, the command
```
./feyntrop < low_level_input.json 2> /dev/null
```
instead.

Changing the number of CPUs used
--------------------------------

By default, **feyntrop** uses the maximal number of available CPUs in the sampling step. This can be changed using the environment variable `OMP_NUM_THREADS`.
For instance, the command
```
OMP_NUM_THREADS=2 ./feyntrop < low_level_input.json 2> /dev/null
```
performs the sampling step for the integration of the input from `low_level_input.json` with only two threads.

A similar option can be used with the python interface.  For instance,
```
OMP_NUM_THREADS=2 python 5L_2pt.py
```
runs the 5 loop 2-point integral example with 2 threads.
