feyntrop
========

Download
--------

 
To download *feyntrop* and the necessary submodules ([Eigen](//eigen.tuxfamily.org/) and [pybind11](//github.com/pybind/pybind11)) use the command

```
git clone --recursive git@github.com:michibo/feyntrop.git
cd feyntrop
```

or 

```
git clone git@github.com:michibo/feyntrop.git
cd feyntrop
git submodule update --init --recursive
```

Compilation
-----------

Make sure to have [python](//python.org) installed.

To compile *feyntrop* call

```
make
```

If the compilation fails, check that the submodules (Eigen and pybind11) are loaded (see instructions above).

Tutorial
--------

To run the tutorial, make sure to have [jupyter-notebook](//jupyter.org/) installed. Then call
```
jupyter notebook tutorial_2L_3pt.ipynb
```
