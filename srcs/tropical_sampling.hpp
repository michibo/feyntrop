/* 

Copyright (c) 2020-2023   Michael Borinsky

*/


// This file implements the heart of the tropical sampling algorithm described in arXiv:2008.12310

#pragma once

#include "graph.hpp"
#include "components.hpp"
#include "random.hpp"

#include "laplacian.hpp"

// cut edge, Jr, r, mm, L
using cut_tuple = tuple< int, double, double, bool, int >;

// Jr, r, mm, L
using J_vector = vector< tuple< double, double, char, char > >;


// Implements the superficial degree of divergence of a subgraph of a Feynman graph
double omega( const graph& g, double D, int L, const edge_subgraph_type& subgraph )
{
    assert( g.is_edge_subgraph(subgraph) );

    double m = 0.0;
    for( int j = 0; j < g._E; j++ )
    {
        if( !subgraph[j] )
            continue;

        pair<int, int> edge;
        double c;
        tie(edge, c) = g._edges[j];

        m += c;
    }

    return m - D/2.0 * L;
}

bool is_mass_momentum_spanning_subgraph(
        const graph& g, 
        const edge_subgraph_type& subgraph, 
        const Eigen::MatrixXd& PGP,
        const Eigen::VectorXd& masses_sqr,
        const vector< Eigen::VectorXd >& Xs,
        const vector< Eigen::VectorXd >& Xinvs,
        const vector<int>& components_map,
        int cV,
        double pmsqr_ref,
        Eigen::MatrixXd& La,
        Eigen::LDLT< Eigen::MatrixXd >& ldlt
        )
{
    assert( g.is_edge_subgraph(subgraph) );
    assert( Xs.size() == Xinvs.size() );
    assert( PGP.cols() == cV-1 && PGP.rows() == cV-1 );
    assert( masses_sqr.size() == g._E );

    La.resize( cV - 1, cV - 1 );

    double sum = 0;
    for( size_t i = 0; i < Xs.size(); i++ )
    {
        get_reduced_contracted_laplacian( La, g, subgraph, cV, components_map, Xinvs[i] );
        ldlt.compute( La );

        double pxi = - ( ldlt.solve( PGP ) ).trace();
        for( int j = 0; j < g._E; j++ )
            if( !subgraph[j] )
                pxi += masses_sqr[j] * Xs[i][j];

        sum += fabs(pxi);
    }

    return sum / pmsqr_ref < 1e-12;
}

bool is_kinematic_generic_3_cut(
        const graph& g, 
        const edge_subgraph_type& subgraph, 
        const Eigen::MatrixXd& PGP,
        const Eigen::VectorXd& masses_sqr,
        const vector<int>& components_map,
        double pmsqr_ref,
        double generic_threshold
        )
{
    assert( PGP.cols() == 2 && PGP.rows() == 2 );

    if( PGP.isZero() )
        return true;

    double a = PGP(0,0);
    double b = PGP(1,1);
    double gamma = PGP(0,1);
    double c = a + b + 2*gamma;

    double A = 0;
    double B = 0;
    double C = 0;

    for( int j = 0; j < g._E; j ++ )
    {
        if( subgraph[j] )
            continue;

        pair<int, int> edge;
        int k,l;
        double c;
        tie(edge, c) = g._edges[j];
        tie(k,l) = edge;

        assert( k != l );
        assert( k > l );

        int kc = components_map[k];
        int lc = components_map[l];

        if( kc == lc )
            continue;

        if( kc == 0 || lc == 0 )
            A += masses_sqr[j];

        if( kc == 1 || lc == 1 )
            B += masses_sqr[j];

        if( kc == 2 || lc == 2 )
            C += masses_sqr[j];
    }

    if( fabs(a - A) / pmsqr_ref > 1e-10 && fabs(b - B) / pmsqr_ref < generic_threshold && fabs(c - C) / pmsqr_ref < generic_threshold )
        return false;

    if( fabs(b - B) / pmsqr_ref > 1e-10 && fabs(a - A) / pmsqr_ref < generic_threshold && fabs(c - C) / pmsqr_ref < generic_threshold )
        return false;

    if( fabs(c - C) / pmsqr_ref > 1e-10 && fabs(a - A) / pmsqr_ref < generic_threshold && fabs(b - B) / pmsqr_ref < generic_threshold )
        return false;

    return true;
}

// checks if Phi/F's Newton polytope is a generalized permutahedron subgraph by subgraph.
// see: Theorem 12.1 of arXiv:1709.07504 on the condition being used here
bool has_subgraph_gp_property(
        const graph& g, 
        const edge_subgraph_type& subgraph, 
        const J_vector& subgraph_table,
        int L,
        bool mm
        )
{
    assert( subgraph.size() >= 2 );

    int z = mm ? L + 1 : L;

    for( int a = 0; a < g._E; a++ )
    {
        if(!subgraph[a])
            continue;

        edge_subgraph_type A = subgraph;
        A.reset(a);

        bool mmAB, mmA, mmB;
        int LAB, LA, LB;
        double Jr, r;
        
        tie( Jr, r, mmA, LA ) = subgraph_table[A.data()];

        int zA = mmA ? LA + 1 : LA;

        for( int b = a+1; b < g._E; b++ )
        {
            if(!subgraph[b])
                continue;

            edge_subgraph_type AB = A; 
            edge_subgraph_type B = subgraph; 
            AB.reset(b);
            B.reset(b);

            tie( Jr, r, mmB, LB ) = subgraph_table[B.data()];
            tie( Jr, r, mmAB, LAB ) = subgraph_table[AB.data()];

            int zB = mmB ? LB + 1 : LB;
            int zAB = mmAB ? LAB + 1 : LAB;

            if( zA + zB > z + zAB )
            {
/*                          // This code provides very detailed information on the violation of the generalized permutahedron property. 
                // It is probably not helpful for people who don't know or care about generalized permutahedra.

                stringstream s;
                s << "Warning: the graph " << g << " Phi/F polynomial's Newton polytope is NOT a generalized permutahedron. With the subgraphs ";
                s << ", A: ";
                for( int j = 0; j < g._E; j++ )
                    if(A[j])
                        s << j << " ";
                s << ", B: ";
                for( int j = 0; j < g._E; j++ )
                    if(B[j])
                        s << j << " ";
                s << ", AnB: ";
                for( int j = 0; j < g._E; j++ )
                    if(AB[j])
                        s << j << " ";
                s << "AuB: ";
                for( int j = 0; j < g._E; j++ )
                    if(subgraph[j])
                        s << j << " ";
                s << ". We get the violation of the property z(A) + z(B) <= z(AnB) + z(AuB).";
                cout << s.str() << endl;
*/
                return false;
            }
        }
    }

    return true;
}

template < class CutFunc >
void visit_cuts( 
        const graph& g, 
        const edge_subgraph_type& subgraph, 
        const J_vector& subgraph_table,
        CutFunc func 
    )
{
    assert( g.is_edge_subgraph(subgraph) );

    for( int j = 0; j < g._E; j++ )
    {
        if( !subgraph[j] )
            continue;

        edge_subgraph_type cutgraph = subgraph;
        cutgraph.reset(j);
       
        double Jr, r;
        bool mm;
        int L;
        tie(Jr, r, mm, L) = subgraph_table[cutgraph.data()];

        if( !func( j, Jr, r, mm, L ) )
            break;
    }
}

double subgraph_Jr_sum( 
        const graph& g, 
        const edge_subgraph_type& subgraph, 
        const J_vector& subgraph_table 
        )
{
    assert( g.is_edge_subgraph(subgraph) );

    double JJr = 0.;
    visit_cuts( g, subgraph, subgraph_table, 
        [&JJr]( int cut_edge, double Jr, double r, bool mm, int L )
        {
            JJr += Jr / r;
            
            return true;
        }
    );

    return JJr;
}

/*

The following function (with the helper functions above) recursively generates a table of the auxillary J_r table as described in Definition 28 and Proposition 29 of arXiv:2008.12310 for the Feynman integral case (see also Section 7.2 of arXiv:2008.12310 )

*/
tuple< J_vector, double, double, bool, bool, bool > generate_subgraph_table( 
        const graph& g, 
        double D, 
        const Eigen::MatrixXd& scalarproducts,
        const Eigen::VectorXd& masses_sqr,
        bool eps_expansion,
        bool bNonEuclidean
        )
{
    assert( D > 0 );
    assert( g._E <= (int)(CHAR_BIT * sizeof(unsigned long long)) );

    if( g._E > (int)(CHAR_BIT * sizeof(unsigned long long) - 1) )
    {
        stringstream s;
        s << "E = " << g._E << " this graph is too large for this implementation!";
        throw domain_error(s.str());
    }

    const double generic_threshold = 1e-5;

    J_vector subgraph_table( 1ULL << g._E );

    vector<int> components_map( g._V );

    edge_subgraph_type cplt_subgraph = g.complete_edge_subgraph();

    int C = components( components_map, g, cplt_subgraph );
    int L = g._E - g._V + C;
    
    if( C != 1 )
    {
        stringstream s;
        s << "Error: the graph " << g << " is not connected.";
        throw domain_error(s.str());
    }

    double W = omega( g, D, L, cplt_subgraph );
    double W_ref = omega( g, -D, L, cplt_subgraph ); // Computes \sum_e \nu_e + D/2 * L
    assert( W_ref > 1e-12 );

    Eigen::MatrixXd PGP;

    edge_subgraph_type subgraph = g.empty_edge_subgraph();
    int cV = components( components_map, g, subgraph );

    assert( g._V == cV );

    get_PGP_matrix( PGP, g, scalarproducts, cV, components_map );
    double pmsqr_ref = PGP.norm() + masses_sqr.norm();

    bool bPhiTrCompute = eps_expansion || fabs(W / W_ref) > 1e-12;

    if( bPhiTrCompute && fabs( pmsqr_ref ) < 1e-12 )
    {
        stringstream s;
        s << "The graph " << g << " has superficial degree of divergence " << W << " and/or we need higher orders in epsilon. So, kinematics are necessary. The incoming momenta/masses appear to be too small to get an (numerically effective) non-zero F/Phi polynomial.";
        throw domain_error(s.str());
    }


    Eigen::MatrixXd La;
    Eigen::LDLT< Eigen::MatrixXd > ldlt;

    int n_tests = 100;
    vector< Eigen::VectorXd > Xs( n_tests );
    vector< Eigen::VectorXd > Xinvs( n_tests );
    
    if( bPhiTrCompute )
    {
        for( int i = 0; i < n_tests; i++ )
        {
            Xs[i] = 2 * Eigen::VectorXd::Ones( g._E ) + Eigen::VectorXd::Random( g._E );
            Xinvs[i] = Xs[i].cwiseInverse();
        }
    }

    bool bGeneric = true;
    bool bPseudoEuclidean = true;
    bool bGPprop = true;

    for( int n = 0; n <= g._E; n++ )
    {
        edge_subgraph_type subgraph = g.empty_edge_subgraph();

        for ( int j = 0; j < n; j++ )
            subgraph.set(j);

        do
        {
            assert( subgraph.count() == n );
            int cV = components( components_map, g, subgraph );

            bool mm = false;
            if( bPhiTrCompute )
            {
                get_PGP_matrix( PGP, g, scalarproducts, cV, components_map );
                assert( (( PGP - PGP.transpose() ) / pmsqr_ref).isZero() );

                if( cV == 2 )
                {
                    assert( PGP.cols() == 1 && PGP.rows() == 1 );

                    if( PGP(0,0) / pmsqr_ref >= 1e-12 ) // We are not in the pseudo Euclidean regime if this is > 0.
                        bPseudoEuclidean = false;
                }

                if( cV == 3 )
                    bGeneric = bGeneric && is_kinematic_generic_3_cut(g, subgraph, PGP, masses_sqr, components_map, pmsqr_ref, generic_threshold );

                mm = is_mass_momentum_spanning_subgraph(g, subgraph, PGP, masses_sqr, Xs, Xinvs, components_map, cV, pmsqr_ref, La, ldlt);
                assert( n != 0 || !mm );
            }
            else if( n == g._E )
                mm = true;

            double Jr, r;

            int L = subgraph.count() - g._V + cV;

            if( n == 0 )
            {
                Jr = 1.0;
                r = 1.0;
            }
            else
            {
                Jr = subgraph_Jr_sum( g, subgraph, subgraph_table );

                r = omega( g, D, L, subgraph );

                if( mm )
                    r -= W;
            }

            if( r/W_ref < 1e-12 && !( n == g._E ) )
            {
                stringstream s;
                s << "Error: the graph " << g << " has a UV or IR subdivergence. \nE.g. the subgraph with edge numbers { ";
                for( int j = 0; j < g._E; j++ )
                    if(subgraph[j])
                        s << j << " ";
                s << "} or, explicitly, with the edges { ";

                for( int j = 0; j < g._E; j++ )
                {
                    if( !subgraph[j] )
                        continue;

                    pair<int, int> edge;
                    int k,l;
                    double c;
                    tie(edge, c) = g._edges[j];
                    tie(k,l) = edge;
                    s << "[" << k << "," << l << "] ";
                }
                s << "}";

                s << " has an (effective) superficial degree of divergence of " << r << ".";
                throw domain_error(s.str());
            }

            // check if Phi/F's Newton polytope is a generalized permutahedron.
            if( bPhiTrCompute && bNonEuclidean && n >= 2 )
                bGPprop = bGPprop && has_subgraph_gp_property(g, subgraph, subgraph_table, L, mm );

            subgraph_table[subgraph.data()] = make_tuple( Jr, r, mm, L );
        }
        while( subgraph.next_permutation() );
    }

    double IGtr = get<0>(subgraph_table[cplt_subgraph.data()]);

    return make_tuple( move(subgraph_table), W, IGtr, bGPprop, bPseudoEuclidean, bGeneric );
}

/*

The remaining functions implement the sampling part of the algorithm described 
in arXiv:2008.12310. I.e. Algorithm 4 of this paper for the Feynman integral 
case. Additionnally, the values of psi^tr and phi^tr are computed.

*/

template< class Generator >
cut_tuple get_random_cut( 
        const graph& g, 
        const edge_subgraph_type& subgraph, 
        double JJr, 
        const J_vector& subgraph_table,
        Generator& gen 
        )
{
    assert( g.is_edge_subgraph(subgraph) );

    double R = true_random::uniform(gen) * JJr;
    double T = .0;

    cut_tuple cut;

    visit_cuts( g, subgraph, subgraph_table, 
        [R,&T,&cut]( int cut_edge, double Jr, double r, bool mm, int L )
        {
            T += Jr / r;

            cut = make_tuple( cut_edge, Jr, r, mm, L );

            if( T > R )
                return false;

            return true;
        }
    );

    assert( T > R );

    return cut;
}

template <class Generator>
pair< double, double > get_random_psi_xi_tropical_sample( 
        Eigen::VectorXd& X, 
        const graph& g, 
        const J_vector& subgraph_table, 
        Generator& gen 
        )
{
    assert( X.size() == g._E );

    edge_subgraph_type subgraph = g.complete_edge_subgraph();

    int L;
    bool MM;
    double JJr, r;

    tie(JJr, r, MM, L) = subgraph_table[subgraph.data()];

    int oL = L;

    double G  = 0.;
    double F  = 1.;
    double XX = 1.;

    for( int j = 0;; j++ )
    {
        int e;
        bool mm;

        tie( e, JJr, r, mm, L ) = get_random_cut( g, subgraph, JJr, subgraph_table, gen );
        assert( subgraph[e] );

        subgraph.reset(e);

        assert( !( mm && !MM ) );

        if( MM && !mm )
            G = XX;

        MM = mm;

        X[e] = XX;

        if( L < oL )
        {
            F *= XX;
            oL = L;
        }

        if( j >= g._E-1 )
            break;

        assert( r > 1e-12 );

        double R = true_random::uniform( gen );
        double Y = pow( 1. - R, 1. / r );

        XX *= Y;
    }

    return make_pair(F, G*F);
}

