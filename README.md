feyntrop
========

**feyntrop** is a computer program to evaluate Feynman integrals. The core C++ integration code, written mainly by [Michael Borinsky](//michaelborinsky.com), is an update of the proof-of-concept implementation [tropical-feynman-quadrature](https://github.com/michibo/tropical-feynman-quadrature), which was published with the paper ['Tropical Monte Carlo quadrature for Feynman integrals'](//arxiv.org/abs/2008.12310). **feyntrop** can be used through a high-level [python](//python.org) interface, written by Henrik Munch.

If **feyntrop** is helpful in your research, please cite,
[M. Borinsky, H. J. Munch, F. Tellander: 'Tropical Feynman integration in the Minkowski regime' arXiv:2302.08955](//arxiv.org/abs/2302.08955) as well as [M. Borinsky: 'Tropical Monte Carlo quadrature for Feynman integrals' arXiv:2008.12310](//arxiv.org/abs/2008.12310).

The implementation internally uses [Eigen](http://eigen.tuxfamily.org), OpenMP and the [xoshiro256+](http://prng.di.unimi.it/) random number generator.
The high-level interface is written in [python](//python.org).

Download
--------

 
To download **feyntrop** and the necessary submodule ([Eigen](//eigen.tuxfamily.org/)) use the command

```
git clone --recursive https://github.com/michibo/feyntrop.git
cd feyntrop
```

or 

```
git clone https://github.com/michibo/feyntrop.git
cd feyntrop
git submodule update --init --recursive
```

Compilation
-----------

To (re)compile **feyntrop** call

```
make clean && make
```

If the compilation fails, check that the submodule (Eigen) is loaded (see instructions above).

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

To run the tutorial in notebook form, you have to have [jupyter-notebook](//jupyter.org/) installed. Run
```
jupyter notebook tutorial_2L_3pt.ipynb
```
in the top directory of this repository to open the tutorial notebook.

Low level interface
-------------------

**feyntrop** can also be used without python. For instance, in a high-performance computing environment. To use this interface, create a file similar to the `low_level_input.json` file in this repository. Here is the content of this file:
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

The content of the json file must be piped into **feyntrop**. For instance, like this:
```
feyntop < low_level_input.json
```
Among some logging information (via stderr), this command produces the output (via stdout) in json format 
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
feyntop < low_level_input.json 2> /dev/null
```
instead.
