/* 

Copyright (c) 2020-2023   Michael Borinsky

*/


#include <tuple>
#include <chrono>

#include <Eigen/Eigenvalues> 

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "graph.hpp"
#include "tropical_sampling.hpp"
#include "feynman_integral.hpp"

int main(int argc, char *argv[])
{
    json data = json::parse(cin);
    //clog << data << endl;

    if( !data.contains("N") || !data.contains("graph") || !data.contains("dimension") || !data.contains("scalarproducts") || !data.contains("masses_sqr") || !data.contains("num_eps_terms") || !data.contains("lambda") )
        throw std::invalid_argument("input json must contain the data fields: N, graph, dimension, scalarproducts, masses_sqr, num_eps_terms, lambda");

    if( !data["N"].is_number() )
        throw std::invalid_argument("input N must be a number");
    if( !data["lambda"].is_number() )
        throw std::invalid_argument("input lambda must be a number");
    if( !data["dimension"].is_number() )
        throw std::invalid_argument("input dimension must be a number");
    if( !data["num_eps_terms"].is_number() )
        throw std::invalid_argument("input num_eps_terms must be a number");

    uint64_t N = data["N"];
    double lambda = data["lambda"];
    double D = data["dimension"];
    int num_eps_terms = data["num_eps_terms"];

    if( !data["graph"].is_array() )
        throw std::invalid_argument("input graph must be a list of the form [[[v1_1,v2_1],nu_1],[[v1_2,v2_2],nu_2],...]");

    graph::edge_vector edges;
    for( auto& edge : data["graph"] )
    {
        if( edge.size() != 2 || edge[0].size() != 2 )
            throw std::invalid_argument("input graph must be a list of the form [[[v1_1,v2_1],nu_1],[[v1_2,v2_2],nu_2],...]");

        if( !edge[0][0].is_number() || !edge[0][1].is_number() || !edge[1].is_number() )
            throw std::invalid_argument("input graph must be a list of the form [[[v1_1,v2_1],nu_1],[[v1_2,v2_2],nu_2],...]");

        int k = edge[0][0];
        int l = edge[0][1];
        double nu = edge[1];

        edges.push_back( make_pair( make_pair(k,l), nu ) );
    }

    graph g( edges );

    if( !data["scalarproducts"].is_array() || data["scalarproducts"].size() != g._V )
        throw std::invalid_argument("input scalarproducts must be a (V x V)-matrix");

    Eigen::MatrixXd scalarproducts(g._V, g._V);

    for( int k = 0; k < g._V; k++ )
    {
        if( !data["scalarproducts"][k].is_array() || data["scalarproducts"][k].size() != g._V )
            throw std::invalid_argument("input scalarproducts must be a (V x V)-matrix");

        for( int l = 0; l < g._V; l++ )
        {
            if( !data["scalarproducts"][k][l].is_number() )
                throw std::invalid_argument("input scalarproducts must be a (V x V)-matrix");
            
            scalarproducts(k,l) = data["scalarproducts"][k][l];
        }
    }

    Eigen::VectorXd masses_sqr( g._E );

    if( !data["masses_sqr"].is_array() || data["masses_sqr"].size() != g._E )
        throw std::invalid_argument("input masses_sqr must be an array of size E");

    for( int e = 0; e < g._E; e++ )
    {
        if( !data["masses_sqr"][e].is_number() )
            throw std::invalid_argument("input masses_sqr must be an array of size E");
            
        masses_sqr[e] = data["masses_sqr"][e];
    }

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
        s << "one mass must be given for each edge.";
        throw domain_error(s.str());
    }

    // clog << "pipj: " << endl;
    // clog << scalarproducts << endl;

    // compute the negative absolute of the scalarproducts matrix
    // i.e. make all eigenvalues negative while keeping Eigen vectors.
    // in the Euclidean case nothing changes as all Eigen values are 
    // negative.
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es( scalarproducts );
    Eigen::MatrixXd scalarproducts_abs = - es.eigenvectors() * es.eigenvalues().cwiseAbs().asDiagonal() * es.eigenvectors().inverse();

    bool bEuclidean = true;
    if( !scalarproducts.isZero() && !scalarproducts.isApprox( scalarproducts_abs ) ) // activate analytic continuation if scalarproducts is not negative-semi-definite
        bEuclidean = false;

    J_vector subgraph_table;

    int W; // superficial degree of divergence w(G)
    double IGtr; // tropicalized Feynman integral I_G^tr
    bool bGeneric;
    bool bPseudoEuclidean;
    bool bGPproperty;

    // time tracking
    chrono::time_point<chrono::system_clock> start, end;
    chrono::duration<double> elapsed_seconds_pre;
    chrono::duration<double> elapsed_seconds_int;


    // Start with Jr subgraph table: (i.e. preprocessing step)
    // clog << "Started calculating Jr-table" << endl;

    start = std::chrono::system_clock::now();
    tie(subgraph_table, W, IGtr, bGPproperty, bPseudoEuclidean, bGeneric) = generate_subgraph_table( g, D, scalarproducts, masses_sqr, num_eps_terms > 1, !bEuclidean );
    end = std::chrono::system_clock::now();

    elapsed_seconds_pre = end-start;

    // clog << "Finished calculating Jr-table in " << elapsed_seconds_pre.count() << " seconds " << endl;
    // clog << "Using " << subgraph_table.size() * sizeof(J_vector::value_type) << " bytes of RAM " << endl;

    // Tropical results:
    // clog << "Superficial degree of divergence: " << W << endl;
    // clog << "I^tr = " << IGtr << endl;

    clog << "(Effective) kinematic regime: ";
    if( bEuclidean )
        clog << "Euclidean";
    else if( bPseudoEuclidean )
        clog << "Pseudo-Euclidean";
    else
        clog << "Minkowski";

    if( bGeneric )
        clog << " (generic).";
    else
        clog << " (exceptional).";
    clog << endl;

    if( bEuclidean )
        assert( bPseudoEuclidean );

    if( bGeneric && !bGPproperty )
    {
        stringstream s;
        s << "The momenta seem generic, but the generalized permutahedron property is not fulfilled. This is likely a bug and the result of the integration shouldn't be trusted. Please contact the developers and send them the following data such that they can fix this bug. Thank you!" << endl;
        s << g << endl;
        s << masses_sqr << endl;
        s << scalarproducts << endl;
        throw domain_error(s.str());
    }

    if( bEuclidean && !bGPproperty )
    {
        stringstream s;
        s << "The momenta seem Euclidean, but the generalized permutahedron property is not fulfilled. This is likely a bug and the result of the integration shouldn't be trusted. Please contact the developers and send them the following data such that they can fix this bug. Thank you!" << endl;
        s << g << endl;
        s << masses_sqr << endl;
        s << scalarproducts << endl;
        throw domain_error(s.str());
    }

    if( !bGeneric && !bEuclidean && W > 0 )
    {   
        //As the integral has positive superficial degree of divergence, we need control over the F polynomials' Newton polytope. 
        clog << "Warning: Kinematics are non-Euclidean and (very) exceptional. Detailed info on the N[F] polytope is needed. The integration might fail or be unstable. Check the result by 1) varying the number of sample points 2) evaluating at multiple kinematic points close to the current one (i.e. by making the kinematics more generic) 3) increasing the spacetime dimension D0." << endl;
        if( bGPproperty )
            clog << "The generalized permutahedron property seems fulfilled. This is a good sign that the integration might succeed anyway." << endl;
        else
            clog << "The generalized permutahedron property does NOT seem fulfilled. Likely, there are singular points in the integration domain that impeed convergence." << endl;
    }
    else
    {
        if( bGPproperty )
        {
            clog << "Generalized permutahedron property seems fulfilled." << endl;
        }
        else
        {
            if( W > 0 )
            {
                stringstream s;
                s << "The momenta seem Euclidean or generic, but the generalized permutahedron property seems not fulfilled. This is likely a bug and the result of the integration shouldn't be trusted. Please contact the developers and send them the following data such that they can fix this bug. Thank you!" << endl;
                s << g << endl;
                s << masses_sqr << endl;
                s << scalarproducts << endl;
                throw domain_error(s.str());
            }
            else
                clog << "Generalized permutahedron property seems NOT fulfilled. However, the integration should work fine as the N[F] polytope is not relevant for the given dimension and edge weights." << endl;
        }
    }

    double deformation_lambda = 0.0;
    if( !bPseudoEuclidean )
    {
        clog << "Analytic continuation: activated. Lambda = " << lambda << "." << endl;

        deformation_lambda = lambda;
    }

    // Initialize multithreading
    int max_threads = omp_get_max_threads();
    omp_set_num_threads(max_threads);

    // Initialize random number generator
    true_random::xoshiro256 gen( 0 ); // Pick your favorite random seed.

    clog << "Started integrating using " << max_threads << " threads and N = " << (double)N << " points." <<  endl;

    start = std::chrono::system_clock::now();
    vector< pair< stats, stats > > res = feynman_integral_estimate( N, g, D, scalarproducts, masses_sqr, num_eps_terms, deformation_lambda, subgraph_table, gen );
    end = std::chrono::system_clock::now();

    elapsed_seconds_int = end-start;

    // Some performance statistics
    clog << "Finished in " << (elapsed_seconds_int + elapsed_seconds_pre).count() << " seconds = " << (elapsed_seconds_int + elapsed_seconds_pre).count()/3600 << " hours." << endl;
    // clog << "Average speed: " << N/elapsed_seconds.count() << " samples / second " << endl;
    // clog << "Relative accuracy: " << res_r.acc()/res_r.avg() << endl;

    // Tropically accelerated Monte Carlo results:
    //clog << "I = (" << res_r.avg() << " + i * " << res_i.avg() << ") +/- (" << res_r.acc() << " + i * " << res_i.acc() << ")" << endl;

    vector< pair< pair< double, double >, pair< double, double > > > res_nums( num_eps_terms );
    for( int j = 0; j < num_eps_terms; j++ )
    {
        res_nums[j] = make_pair( make_pair( get<0>(res[j]).avg(), get<0>(res[j]).acc()), 
                                 make_pair( get<1>(res[j]).avg(), get<1>(res[j]).acc()) );
    }

    json output;

    output["integral"] = res_nums;
    output["IGtr"] = IGtr;
    output["seconds preprocessing"] = elapsed_seconds_pre.count();
    output["seconds sampling"] = elapsed_seconds_int.count();

    cout << output.dump() << endl;

    return 0;
}

