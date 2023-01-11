/* 

Copyright (c) 2020-2023   Michael Borinsky

*/

#include <complex>

#include "stats.hpp"

#include "random.hpp"
#include "tropical_sampling.hpp"
#include "symanzik_polynomials.hpp"

void infinite_value_warning( const graph& g, complex<double> R, complex<double> cPsi, complex<double> cPPhi, const Eigen::VectorXd& X, const Eigen::VectorXcd& cX )
{
    #pragma omp critical
    {
        cout << "Warning: Sampled value " << R << " - floating point accuracy or numerical stabilty seem to be insufficient - dropping this point" << endl;
        cout << "If this happens often, the result will not be reliable!" << endl;

        for( int j = 0; j < g._E; j++ )
            cout << "x_" << j << " = " << X[j] << "; ";
        cout << endl;
        for( int j = 0; j < g._E; j++ )
            cout << "cX_" << j << " = " << cX[j] << "; ";
        cout << endl;

        cout << "cPsi = " << cPsi << " ; cPPhi = " << cPPhi << endl;
        cout << "-" << endl;
    }
}

template< class Generator >
vector< pair< stats, stats > > feynman_integral_estimate( 
        uint64_t N, 
        const graph& g, 
        double D, 
        const Eigen::MatrixXd& scalarproducts, 
        const Eigen::MatrixXd& scalarproducts_abs, 
        const Eigen::VectorXd& masses_sqr, 
        int num_eps_terms, 
        double deformation_lambda, 
        const J_vector& subgraph_table, 
        Generator& gen
        )
{
    assert( scalarproducts.cols() == g._V || scalarproducts.rows() == g._V );
    assert( scalarproducts.isApprox( scalarproducts.transpose() ) );
    assert( masses_sqr.size() == g._E );
    assert( subgraph_table.size() == ( 1ULL << g._E ) );
    assert( D > 1e-12 );

    using namespace Eigen;

    edge_subgraph_type subgraph = g.complete_edge_subgraph();

    double IGtr;
    bool MM;
    int L;
    double r;
    tie(IGtr, r, MM, L) = subgraph_table[subgraph.data()];
    assert(MM);

    double W = omega( g, D, L, subgraph );
    double W_ref = omega( g, -D, L, subgraph ); // Computes \sum_e \nu_e + D/2 * L
    assert( W_ref > 1e-12 );

    bool W_zero = fabs( W/W_ref ) < 1e-12;

    edge_subgraph_type contracted_subgraph = g.empty_edge_subgraph();
    vector<int> contracted_subgraph_components_map( g._V );
    int C = components( contracted_subgraph_components_map, g, contracted_subgraph );

    assert( C == g._V );

    VectorXd nu( g._E );
    for( int j = 0; j < g._E; j++ )
    {
        pair<int, int> edge;
        double c;
        tie(edge, c) = g._edges[j];

        nu[j] = c;
    }

    MatrixXd PGP, PGP_abs;
    get_PGP_matrix( PGP, g, scalarproducts, C, contracted_subgraph_components_map );
    get_PGP_matrix( PGP_abs, g, scalarproducts_abs, C, contracted_subgraph_components_map );

    MatrixXcd cPGP = PGP;

    vector< true_random::xoshiro256 > generators; 

    int max_threads = omp_get_max_threads();

    for( int i = 0; i < max_threads; i++ )
    {
        generators.emplace_back( gen );
        gen.jump();
    }

    vector< pair< vector<stats>, vector<stats> > > mcs( num_eps_terms );

    for( int j = 0; j < num_eps_terms; j++ )
        mcs[j] = make_pair( vector<stats>( max_threads ), vector<stats>( max_threads ) );

    MatrixXd La( g._V-1, g._V-1 );
    MatrixXd LaInv( g._V-1, g._V-1 );
    LDLT< MatrixXd > ldlt( g._V-1 );

    VectorXd X( g._E );
    VectorXd Xinv( g._E );
    VectorXd XinvSqr( g._E );

    VectorXd d_pphi( g._E );
    MatrixXd dd_pphi( g._E, g._E );

    MatrixXd LPGPLT( g._V-1, g._V-1 );

    VectorXcd cX( g._E );
    VectorXcd cXinv( g._E );
    MatrixXcd cLa( g._V-1, g._V-1 );

    FullPivLU<MatrixXcd> luLa( g._V-1, g._V-1 );

    MatrixXcd cJacobian( g._E, g._E );
    //FullPivLU<MatrixXcd> luJac( g._E, g._E );
    PartialPivLU<MatrixXcd> luJac( g._E );

    double def_ref = 0.;
    if( deformation_lambda != 0 )
    {
        Xinv = VectorXd::Ones( g._E );
        XinvSqr = VectorXd::Ones( g._E );

        get_reduced_weighted_laplacian( La, g, Xinv );
        ldlt.compute( La );
        LaInv = ldlt.solve( MatrixXd::Identity( g._V-1, g._V-1 ) );
        LPGPLT = LaInv * PGP_abs * LaInv.transpose();

        eval_pphi_polynomial_derivatives( d_pphi, dd_pphi, g, Xinv, XinvSqr, LaInv, LPGPLT );

        def_ref = -M_PI/2. * deformation_lambda / ( d_pphi.sum() + masses_sqr.sum() );
    }
    else
    {
        for( int t = 0; t < max_threads; t++ )
        {
            for( int j = 0; j < num_eps_terms; j++ )
            {
                get<1>(mcs[j])[t].update(0);
                get<1>(mcs[j])[t].update(0);
            }
        }
    }

    #pragma omp parallel for default(none) \
        shared(cout,subgraph_table,mcs,generators,N,g,D,L,W,W_zero,PGP,cPGP,masses_sqr,nu,num_eps_terms,def_ref,IGtr) \
        firstprivate(La,LaInv,ldlt,X,Xinv,XinvSqr,d_pphi,dd_pphi,LPGPLT,cX,cXinv,cLa,luLa,cJacobian,luJac) \
        schedule(dynamic, 10000)
    for( uint64_t i = 0; i < N; i++ )
    {
        int t = omp_get_thread_num();

        // generate the random sample from the tropicalized measure of the Feynman integral and the values of psi_tr and xi_tr
        double psi_tr, xi_tr;
        tie(psi_tr, xi_tr) = get_random_psi_xi_tropical_sample( X, g, subgraph_table, generators[t] );

        /// !

        // it is good for numerical stability to normalize the X variables 
        // such that psi_tr^(-D/2) * (psi_tr/xi_tr)^W = 1

        double tgt = pow(psi_tr, -D/2.) * pow(psi_tr/xi_tr, W);
        double scale = pow( tgt, 1. / ( D/2. * L + W ) );

        X *= scale;

        psi_tr = 1.0;
        xi_tr = 1.0;

        /// ! remove ! block if this is not needed
        /// !!! New: If this block is removed also later code needs to be changed...

        Xinv = X.cwiseInverse(); // coefficient wise inverse of X variables.

        get_reduced_weighted_laplacian( La, g, Xinv );

        ldlt.compute( La );

        if( def_ref == 0. || ( W_zero && num_eps_terms == 1 ) )
        {
            double detCholeskyL = ldlt.matrixL().determinant();
            double detCholeskyD = ldlt.vectorD().prod();

            double psi = detCholeskyL * detCholeskyL * detCholeskyD * X.prod();
            
            //double R = IGtr * pow( psi_tr/psi, D/2. );
            double R = IGtr * pow( psi, -D/2. );

            if( !isfinite( R ) )
            {
                infinite_value_warning( g, R, psi, 0., X, cX );
                continue;
            }

            // don't need to compute second Symanzik if W = 0!
            if( W_zero && num_eps_terms == 1 )
            {
                get<0>(mcs[0])[t].update( R );
                continue;
            }

            double pphi = - ( ldlt.solve( PGP ) ).trace();
            double pxi  = pphi + masses_sqr.dot(X);

            //double xi  = psi * ( pphi + masses_sqr.dot(X) );

            //R *= pow( (psi / xi ) / (psi_tr / xi_tr), W );
            R *= pow( pxi, -W );

            if( !isfinite( R ) )
            {
                infinite_value_warning( g, R, psi, pphi, X, cX );
                continue;
            }

            get<0>(mcs[0])[t].update( R );

            if( num_eps_terms > 1 )
            {
                double eps_fac = log(psi) - L * log(pxi);
                for( int j = 1; j < num_eps_terms; j++ )
                {
                    R *= eps_fac / (double)j;
                    get<0>(mcs[j])[t].update( R );
                }
            }

            continue;
        }

        // Analytic continuation starts here

        LaInv = ldlt.solve( MatrixXd::Identity( g._V-1, g._V-1 ) );
        LPGPLT = LaInv * PGP * LaInv.transpose();
        XinvSqr = Xinv.array().square(); // coefficient wise square

        eval_pphi_polynomial_derivatives( d_pphi, dd_pphi, g, Xinv, XinvSqr, LaInv, LPGPLT );

        // pphi derivative gets modified to hold the deformation phase vector
        d_pphi += masses_sqr;
        d_pphi *= def_ref;

        for( int j = 0; j < g._E; j++ )
            cX[j] = X[j] * complex<double>( cos( d_pphi[j] ), sin( d_pphi[j] ) );

        cXinv = cX.cwiseInverse();

        get_reduced_weighted_laplacian( cLa, g, cXinv );

        // The Laplacian is only half-filled. The upper half is filled in this loop
        for( int k = 0; k < g._V-1; k++ )
            for( int l = k+1; l < g._V-1; l++ )
                cLa( k, l ) = cLa( l, k );
    
        luLa.compute( cLa );

        complex<double> cPsi = cX.prod() * luLa.determinant();

        //complex<double> cR = IGtr * pow( psi_tr/cPsi, D/2. );
        complex<double> cR = IGtr * pow( cPsi, -D/2. );

        cJacobian = complex<double>( 0., def_ref ) * dd_pphi * X.asDiagonal() + MatrixXd::Identity( g._E, g._E );

        luJac.compute( cJacobian );

        complex<double> det_jac = luJac.determinant();

        double x_nu_exp = nu.dot(d_pphi);

        cR *= det_jac * complex<double>( cos( x_nu_exp ), sin( x_nu_exp ) );

        complex<double> cPPhi = - ( luLa.solve( cPGP ) ).trace();
        complex<double> cPXi = cPPhi + masses_sqr.dot(cX);

        //complex<double> cXi = cPsi * ( cPPhi + masses_sqr.dot(cX) );
        //cR *= pow( (cPsi / cXi ) / (psi_tr / xi_tr), W );
        cR *= pow( cPXi, -W );

        if( !isfinite( real(cR) ) || !isfinite( imag(cR) ) )
        {
            infinite_value_warning( g, cR, cPsi, cPPhi, X, cX );
            continue;
        }

        get<0>(mcs[0])[t].update( real(cR) );
        get<1>(mcs[0])[t].update( imag(cR) );

        if( num_eps_terms > 1 )
        {
            complex<double> eps_fac = log(cPsi) - (double)L * log(cPXi);
            for( int j = 1; j < num_eps_terms; j++ )
            {
                cR *= eps_fac / (double)j;
                get<0>(mcs[j])[t].update( real(cR) );
                get<1>(mcs[j])[t].update( imag(cR) );
            }
        }
    }

    vector< pair< stats, stats > > result( num_eps_terms );
    for( int j = 0; j < num_eps_terms; j++ )
        result[j] = make_pair( merge_stats_vector( get<0>(mcs[j]) ), merge_stats_vector( get<1>(mcs[j]) ) );

    return result;
}

