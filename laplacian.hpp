/* 

Copyright (c) 2020-2023   Michael Borinsky

*/

#pragma once

#include <vector>
#include <Eigen/Dense>

#include "graph.hpp"

template< class Field>
void get_reduced_weighted_laplacian( 
        Eigen::Matrix<Field, Eigen::Dynamic, Eigen::Dynamic>& La,
        const graph& g, 
        const Eigen::Matrix<Field, Eigen::Dynamic, 1>& Xinv
        )
{
    assert( Xinv.size() == g._E );

    int E = g._E;
    int V = g._V;
    
    assert( La.rows() == V-1 && La.cols() == V-1 );

    La.setZero( V-1, V-1 );

    for ( int j = 0; j < E; j++ )
    {
        pair<int, int> edge;
        int k,l;
        double c;
        tie(edge, c) = g._edges[j];
        tie(k,l) = edge;

        assert( k != l );
        assert( k > l );

        Field x = Xinv[j];

        La(l,l) +=  x;

        if( k == V-1 ) // delete last row
            continue;

        La(k,k) +=  x;
        La(k,l) += -x;
        // only lower triangular part of the laplacian matters for ldlt module. 
        // => Fill only lower triangular part.
        // For the use with a complex Laplacian, we need to fill also the other part!
    }
}

void get_PGP_matrix( 
        Eigen::MatrixXd& PGP,
        const graph& g, 
        const Eigen::MatrixXd& scalarproducts,
        int num_components,
        const vector<int>& contracted_subgraph_components_map
        )
{
    assert( contracted_subgraph_components_map.size() == g._V );

    assert( scalarproducts.cols() == g._V && scalarproducts.rows() == g._V );
    assert( scalarproducts.isApprox( scalarproducts.transpose() ) );

    int V = num_components;
    PGP.setZero( V-1, V-1 );
    
    for( int i = 0; i < g._V; i++ )
    {
        for( int j = 0; j < g._V; j++ )
        {
            int v = contracted_subgraph_components_map[i];
            int w = contracted_subgraph_components_map[j];

            if( v == V-1 || w == V-1 )
                continue;

            PGP(v,w) += scalarproducts(i,j);
        }
    }
}

// computes Phi^0/Psi and its x_e derivatives
void eval_pphi_polynomial_derivatives(
        Eigen::VectorXd& d_pphi,
        Eigen::MatrixXd& dd_pphi,
        const graph& g, 
        const Eigen::VectorXd& Xinv,
        const Eigen::VectorXd& XinvSqr,
        const Eigen::MatrixXd& LaInv,
        const Eigen::MatrixXd& LPGPLT
        )
{
    assert( Xinv.size() == g._E );
    assert( XinvSqr.size() == g._E );
    assert( d_pphi.size() == g._E );
    assert( dd_pphi.cols() == g._E );
    assert( dd_pphi.rows() == g._E );
    assert( LaInv.rows() == g._V-1 );
    assert( LaInv.cols() == g._V-1 );
    assert( LPGPLT.rows() == g._V-1 );
    assert( LPGPLT.cols() == g._V-1 );

    int V = g._V;
    int E = g._E;

    for( int a = 0; a < E; a++ )
    {
        pair<int, int> edge_a;
        int ka,la,ca;
        tie(edge_a, ca) = g._edges[a];
        tie(ka,la) = edge_a;
        assert( ka > la );

        double s, t;
            
        if( ka == V-1 )
        {
            s = LPGPLT(la,la);
            t = LaInv(la,la);
        }
        else
        {
            s = LPGPLT(ka,ka) + LPGPLT(la,la) - 2. * LPGPLT(ka,la);
            t = LaInv(ka,ka) + LaInv(la,la) - 2. * LaInv(ka,la);
        }

        double xa = XinvSqr[a];
        double xas = xa * s;

        d_pphi(a) = - xas;
        dd_pphi(a,a) = 2. * xas * ( Xinv[a] - xa * t );

        for( int b = a + 1; b < E; b++ )
        {
            pair<int, int> edge_b;
            int kb,lb,cb;
            tie(edge_b, cb) = g._edges[b];
            tie(kb,lb) = edge_b;
            assert( kb > lb );

            double s, t;

            if( kb == V-1 )
            {
                if( ka == V-1 )
                {
                    s = LPGPLT(la,lb);
                    t = LaInv(la,lb);
                }
                else
                {
                    s = LPGPLT(la,lb) - LPGPLT(ka,lb);
                    t = LaInv(la,lb) - LaInv(ka,lb);
                }
            }
            else 
            {
                if( ka == V-1 )
                {
                    s = LPGPLT(la,lb) - LPGPLT(la,kb);
                    t = LaInv(la,lb) - LaInv(la,kb);
                }
                else
                {
                    s = LPGPLT(ka,kb) + LPGPLT(la,lb) - LPGPLT(ka,lb) - LPGPLT(la,kb);
                    t = LaInv(ka,kb) + LaInv(la,lb) - LaInv(ka,lb) - LaInv(la,kb);
                }
            }
            
            double xb = XinvSqr[b];

            dd_pphi(a,b) = dd_pphi(b,a) = - 2. * xa * xb * s * t;
        }
    }
}

