A tropical Monte Carlo integrator for parametric Feynman integrals
==================================================================

This repository is a proof-of-concept implementation of an algorithm to estimate Feynman integrals using tropical sampling. The details of this algorithm and its general form for arbitrary generalized permutahedra type integrals is described in the paper ['Tropical Monte Carlo quadrature for Feynman integrals'](//arxiv.org/abs/2008.12310). 

The implementation uses the [xoshiro256+](http://prng.di.unimi.it/) random number generator.

Prerequisites
-------------

The following programs/packages need to be available to run the code:

- gcc
- OpenMP and 
- the [Eigen](http://eigen.tuxfamily.org) linear algebra C++ library

Compiling and running
---------------------

Run `> make` to compile the code. If the compilation was successful, you can run the example programs via 
```
> ./simple_example
``` 

or 

```
> ./advanced_example
```


Example
-------

The algorithm is implemented in the following header files that have to be included:

```C++
#include "graph.hpp"
#include "subgraph_table.hpp"
#include "feynman_integral.hpp"
```

Here is an example code of the integration of the wheel with three spokes graph period integral:

```C++
// We are going to integrate the 
// wheel with 3 spokes graph 
graph g( 
  graph::edge_vector{ 
    graph::edge_tuple{ {0, 1}, 1 }, 
    graph::edge_tuple{ {0, 2}, 1 }, 
    graph::edge_tuple{ {0, 3}, 1 }, 
    graph::edge_tuple{ {1, 2}, 1 }, 
    graph::edge_tuple{ {2, 3}, 1 }, 
    graph::edge_tuple{ {3, 1}, 1 }, 
  } 
);

// Notation for edges:
//
// the object graph::edge_tuple{ {v1, v2}, w }
// represents an edge from vertex v1 to vertex v2
// with edge weight w


vector< Eigen::VectorXd > momenta{
    Eigen::Vector4d{   3.0, - 1.0, - 1.0, - 1.0 },
    Eigen::Vector4d{ - 1.0,   3.0, - 1.0, - 1.0 }, 
    Eigen::Vector4d{ - 1.0, - 1.0,   3.0, - 1.0 }, 
    Eigen::Vector4d{ - 1.0, - 1.0, - 1.0,   4.0 }
  };

// one incoming momentum for each vertex of the graph 
// is needed. An incoming momentum can be 0!

// The momenta must be Euclidean! Minkowski vectors are 
// not implemented (yet)


vector<double> masses_sqr{ 
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0 
  };

// one squared mass for each edge of the graph is needed.


// We consider a D = 4 integral:
constexpr int D = 4;


double mm_eps = 1e-10; 
// mm_eps must be almost zero
// it is used to decide if a subgraph is 
// mass-momentum spanning.


// This table stores the probabilities for each maximal 
// cone of the braid arrangement fan:
J_vector subgraph_table;


int W; // superficial degree of divergence w(G)
double IGtr; // tropicalized Feynman integral I_G^tr


// Compute the Jr subgraph table: (i.e. perform the preprocessing step)
tie(subgraph_table, W, IGtr) = generate_subgraph_table( g, D, momenta, masses_sqr, mm_eps );


// 'Tropical' results:
cout << "Superficial degree of divergence: " << W << endl;
cout << "I^tr = " << IGtr << endl;


// Initialize random number generator
true_random::xoshiro256 gen( 0 );


// Number of points to be sampled:
constexpr uint64_t N = 1000000ULL;


// Perform the actual Monte Carlo integration:
stats res = feynman_integral_estimate( N, g, D, momenta, masses_sqr, subgraph_table, gen );


// the res object stores the result:
//
// res.avg() gives the estimate,
// res.acc() the estimated accuracy
// res.var() the estimated sample variance

// Monte Carlo results:
cout << "I = " << res.avg() << " +/- " << res.acc() << endl;

cout << "Relative accuracy: " << res.acc()/res.avg() << endl;
```

This example is implemented in `simple_example.cpp` and can be executed by running the `> ./simple_example` program, which is compiled using `> make`.

More advanced examples with some benchmark code, integration of massive Feynman graphs and graphs with more loops can be found in `advanced_example.cpp`, which can be run using `> ./advanced_example`. 

