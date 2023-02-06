/* 

Copyright (c) 2020-2023   Michael Borinsky

*/


#include <tuple>
#include <chrono>

#include <Eigen/Eigenvalues> 

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>
#include <pybind11/eigen.h>

#include "graph.hpp"
#include "tropical_sampling.hpp"
#include "feynman_integral.hpp"

namespace py = pybind11;

pair< vector< pair< pair< double, double >, pair< double, double > > >, double > integrate_graph( const graph& g, int D, const Eigen::MatrixXd& scalarproducts, const Eigen::VectorXd& masses_sqr, int num_eps_terms, double lambda, uint64_t N )
{
    if( scalarproducts.cols() != g._V || scalarproducts.rows() != g._V || !scalarproducts.isApprox( scalarproducts.transpose() ) )
    {
        stringstream s;
        s << "integrate_graph: The matrix \n" << scalarproducts << "\nof invariants was given as input and it is not symmetric or not a VxV matrix. One row and column must be specified for each vertex - representing the momentum flowing into this vertex.";
        throw domain_error(s.str());
    }

    if( !scalarproducts.colwise().sum().isZero() )
    {
        stringstream s;
        s << "Scalar product (P) matrix has to fulfill momentum conservation. That means all row and column sums need to be 0. The matrix \n";
        s << scalarproducts << "\n";
        s << "does not fulfill this requirement.";
        throw domain_error(s.str());
    }

    if( masses_sqr.size() != g._E )
    {
        stringstream s;
        s << "integrate_graph: one mass must be given for each edge.";
        throw domain_error(s.str());
    }

    // cout << "pipj: " << endl;
    // cout << scalarproducts << endl;

    // compute the negative absolute of the scalarproducts matrix
    // i.e. make all eigenvalues negative while keeping Eigen vectors.
    // in the Euclidean case nothing changes as all Eigen values are 
    // negative.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es( scalarproducts );
    Eigen::MatrixXd scalarproducts_abs = - es.eigenvectors() * es.eigenvalues().cwiseAbs().asDiagonal() * es.eigenvectors().inverse();
    

    bool bEuclidean = true;
    if( !scalarproducts.isApprox( scalarproducts_abs ) ) // activate analytic continuation if scalarproducts is not negative-semi-definite
        bEuclidean = false;

    J_vector subgraph_table;

    int W; // superficial degree of divergence w(G)
    double IGtr; // tropicalized Feynman integral I_G^tr
    bool bPseudoEuclidean;
    bool bGPproperty;

    // time tracking
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed_seconds;


    // Start with Jr subgraph table: (i.e. preprocessing step)
    // cout << "Started calculating Jr-table" << endl;

    start = std::chrono::system_clock::now();
    tie(subgraph_table, W, IGtr, bGPproperty, bPseudoEuclidean) = generate_subgraph_table( g, D, scalarproducts, masses_sqr, num_eps_terms > 1, !bEuclidean );
    end = std::chrono::system_clock::now();

    elapsed_seconds = end-start;

    // cout << "Finished calculating Jr-table in " << elapsed_seconds.count() << " seconds " << endl;
    // cout << "Using " << subgraph_table.size() * sizeof(J_vector::value_type) << " bytes of RAM " << endl;

    // Tropical results:
    // cout << "Superficial degree of divergence: " << W << endl;
    // cout << "I^tr = " << IGtr << endl;

    if(bGPproperty)
        cout << "Generalized permutahedron property: fulfilled." << endl;
    else
        cout << "Generalized permutahedron property: NOT fulfilled. Integration may fail. Check the result by varying the number of sample points and by evaluating at multiple kinematic points close to the current one." << endl;

    cout << "(Effective) kinematic regime: ";
    if( bEuclidean )
        cout << "Euclidean";
    else if( bPseudoEuclidean )
        cout << "Pseudo-Euclidean";
    else
        cout << "Minkowski";
    cout << endl;

    if( bEuclidean )
        assert( bPseudoEuclidean );

    double deformation_lambda = 0.0;
    if( !bPseudoEuclidean )
    {
        cout << "Analytic continuation: activated. Lambda = " << lambda << endl;

        deformation_lambda = lambda;
    }

    // Initialize multithreading
    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);

    // Initialize random number generator
    true_random::xoshiro256 gen( 0 ); // Pick your favorite random seed.

    cout << "Started integrating using " << max_threads << " threads and N = " << (double)N << " points." <<  endl;

    start = std::chrono::system_clock::now();
    vector< pair< stats, stats > > res = feynman_integral_estimate( N, g, D, scalarproducts, masses_sqr, num_eps_terms, deformation_lambda, subgraph_table, gen );
    end = std::chrono::system_clock::now();

    elapsed_seconds = end-start;

    // Some performance statistics
    cout << "Finished in " << elapsed_seconds.count() << " seconds = " << elapsed_seconds.count()/3600 << " hours." << endl;
    // cout << "Average speed: " << N/elapsed_seconds.count() << " samples / second " << endl;
    // cout << "Relative accuracy: " << res_r.acc()/res_r.avg() << endl;

    // Tropically accelerated Monte Carlo results:
    //cout << "I = (" << res_r.avg() << " + i * " << res_i.avg() << ") +/- (" << res_r.acc() << " + i * " << res_i.acc() << ")" << endl;

    vector< pair< pair< double, double >, pair< double, double > > > res_nums( num_eps_terms );
    for( int j = 0; j < num_eps_terms; j++ )
    {
        res_nums[j] = make_pair( make_pair( get<0>(res[j]).avg(), get<0>(res[j]).acc()), 
                                 make_pair( get<1>(res[j]).avg(), get<1>(res[j]).acc()) );
    }

    return make_pair(res_nums, IGtr);
}

PYBIND11_MODULE(feyntrop, m) {
    m.doc() = "feyntrop plugin";

    py::class_<graph>(m, "graph")
        .def(py::init<graph::edge_vector>())
        .def_readwrite("E", &graph::_E)
        .def_readwrite("V", &graph::_V)
        .def_readwrite("edges", &graph::_edges)
        .def("__repr__",
            [](const graph& g) {
                stringstream s;
                s << g;
                return s.str();
            });

    m.def("integrate_graph", &integrate_graph,
            py::call_guard<py::scoped_ostream_redirect,
                           py::scoped_estream_redirect>()
        );
}

