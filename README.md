feyntrop
========

**feyntrop** is a computer program to evaluate Feynman integrals. The core C++ integration code, written mainly by [Michael Borinsky](//michaelborinsky.com), is an update of the proof-of-concept implementation [tropical-feynman-quadrature](https://github.com/michibo/tropical-feynman-quadrature), which was published with the paper ['Tropical Monte Carlo quadrature for Feynman integrals'](//arxiv.org/abs/2008.12310). **feyntrop** can be used through a [python](//python.org) interface, written by Henrik Munch.

If **feyntrop** is helpful in your research, please cite,
[M. Borinsky, H. J. Munch, F. Tellander: 'Tropical Feynman integration in the Minkowski regime' arXiv:2302.08955](//arxiv.org/abs/2302.08955) as well as [M. Borinsky: 'Tropical Monte Carlo quadrature for Feynman integrals' arXiv:2008.12310](//arxiv.org/abs/2008.12310).

The implementation internally uses [Eigen](http://eigen.tuxfamily.org), [python](//python.org), [pybind11](//github.com/pybind/pybind11)), OpenMP and the [xoshiro256+](http://prng.di.unimi.it/) random number generator.

Download
--------

 
To download **feyntrop** and the necessary submodules ([Eigen](//eigen.tuxfamily.org/) and [pybind11](//github.com/pybind/pybind11)) use the command

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

Make sure to have [python](//python.org) installed.

To compile **feyntrop** call

```
make clean && make
```

If the compilation fails, check that the submodules (Eigen and pybind11) are loaded (see instructions above).

Tests
-----

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
