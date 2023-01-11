/* 

Copyright (c) 2020-2023   Michael Borinsky

*/


#include <array>
#include <cmath>

constexpr int N_SUM_LIMBS = 64;

class fp_summer
{
public:
    fp_summer()
    {
        for( double& s : _s )
            s = .0;
    }

    void add( double v )
    {
        if( v == .0 )
            return;

        bool bAdded = false;
        for( double& s : _s )
        {
            if( s == .0 )
            {
                s = v;

                bAdded = true;
                break;
            }

            v += s;
            s = .0;
        }

        if(!bAdded)
            throw runtime_error("ERROR: More limbs needed for fp_summer");
    }

    double value() const
    {
        double val = .0;

        for( double s : _s )
            val += s;

        return val;
    }

    array< double, N_SUM_LIMBS > _s;
};

fp_summer merge_summers( const fp_summer& entry1, const fp_summer& entry2 )
{
    assert( entry1._s[N_SUM_LIMBS-1] == .0 && entry2._s[N_SUM_LIMBS-1] == .0 );

    fp_summer m;

    for( int i = 0; i < N_SUM_LIMBS-1; i++ )
        m._s[i+1] = entry1._s[i] + entry2._s[i];

    m._s[0] = .0;

    return m;
}

class stats
{
public:
    stats()
    { }

    void update( double val )
    {
        N++;
        S1.add( val );
        S2.add( val*val );
    }

    double avg() const
    {
        return S1.value() / N;
    }

    double var() const
    {
        return ((double)N / ((double)N - 1)) * ( S2.value() / N - pow( avg(), 2 ) );
    }

    double std_def() const
    {
        return sqrt(var());
    }

    double acc() const
    {
        return sqrt(var()/N);
    }

public:
    uint64_t N = 0;

    fp_summer S1, S2;
};

stats merge_stats( const stats& mc1, const stats& mc2 )
{
    stats mc;
    mc.N = mc1.N + mc2.N;

    mc.S1 = merge_summers( mc1.S1, mc2.S1 );
    mc.S2 = merge_summers( mc1.S2, mc2.S2 );

    return mc;
}

stats merge_stats_vector( vector<stats> mcs )
{
    for( size_t e = 0; (1ULL << e) < mcs.size(); e++ )
    {
        for( size_t i = 0; i + (1ULL << e) < mcs.size(); i += 1ULL << (e+1) )
        {
            mcs[i] = merge_stats( mcs[i], mcs[i + (1 << e)] );
        }
    }

    return mcs[0];
}
