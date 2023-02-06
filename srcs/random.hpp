
// Adapted by Michael Borinsky 2020-2023

// uses algorithms from arXiv:1805.01407
// by David Blackman and Sebastiano Vigna
// see http://prng.di.unimi.it/

#pragma once

#include <cmath>
#include <cstdint>
#include <climits>

namespace true_random
{
    inline double uint64_to_double_01( uint64_t x )
    {
        return (x >> 11) * (1. / (UINT64_C(1) << 53));
    }

    // generates with measure: dmu = dx
    // with x in [0,1)
    template<class Generator>
    inline double uniform( Generator& gen )
    {
        return uint64_to_double_01(gen.next());
    }

    static inline uint64_t rotl(const uint64_t x, int k) {
        return (x << k) | (x >> (64 - k));
    }

    class splitmix
    {
    public:
        splitmix( uint64_t seed ) : x( seed )
        { }
       
        /*  Written in 2015 by Sebastiano Vigna (vigna@acm.org)

        To the extent possible under law, the author has dedicated all copyright
        and related and neighboring rights to this software to the public domain
        worldwide. This software is distributed without any warranty.

        See <http://creativecommons.org/publicdomain/zero/1.0/>. */

        /* This is a fixed-increment version of Java 8's SplittableRandom generator
           See http://dx.doi.org/10.1145/2714064.2660195 and 
           http://docs.oracle.com/javase/8/docs/api/java/util/SplittableRandom.html

           It is a very fast generator passing BigCrush, and it can be useful if
           for some reason you absolutely want 64 bits of state; otherwise, we
           rather suggest to use a xoroshiro128+ (for moderately parallel
           computations) or xorshift1024* (for massively parallel computations)
           generator. */

    private:
        uint64_t x; /* The state can be seeded with any value. */

    public:
        uint64_t next() {
            uint64_t z = (x += 0x9e3779b97f4a7c15);
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
            z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
            return z ^ (z >> 31);
        }
    };

    class xoshiro256
    {
    public:

        using result_type = uint64_t;

        static result_type min()
        {
            return 0;
        }

        static result_type max()
        {
            return ULLONG_MAX;
        }

    public:
        xoshiro256( uint64_t seed )
        {
            splitmix m(seed);

            for( int i = 0; i < 4; i++ )
                s[i] = m.next();
        }

        /*  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

        To the extent possible under law, the author has dedicated all copyright
        and related and neighboring rights to this software to the public domain
        worldwide. This software is distributed without any warranty.

        See <http://creativecommons.org/publicdomain/zero/1.0/>. */
        /* This is xoshiro256+ 1.0, our best and fastest generator for floating-point
           numbers. We suggest to use its upper bits for floating-point
           generation, as it is slightly faster than xoshiro256**. It passes all
           tests we are aware of except for the lowest three bits, which might
           fail linearity tests (and just those), so if low linear complexity is
           not considered an issue (as it is usually the case) it can be used to
           generate 64-bit outputs, too.

           We suggest to use a sign test to extract a random Boolean value, and
           right shifts to extract subsets of bits.

           The state must be seeded so that it is not everywhere zero. If you have
           a 64-bit seed, we suggest to seed a splitmix64 generator and use its
           output to fill s. */

    private:
        uint64_t s[4];

    public:
        result_type operator()()
        {
            return next();
        }

        uint64_t next(void) {
            const uint64_t result_plus = s[0] + s[3];

            const uint64_t t = s[1] << 17;

            s[2] ^= s[0];
            s[3] ^= s[1];
            s[1] ^= s[2];
            s[0] ^= s[3];

            s[2] ^= t;

            s[3] = rotl(s[3], 45);

            return result_plus;
        }


        /* This is the jump function for the generator. It is equivalent
           to 2^128 calls to next(); it can be used to generate 2^128
           non-overlapping subsequences for parallel computations. */

        void jump(void) {
            static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

            uint64_t s0 = 0;
            uint64_t s1 = 0;
            uint64_t s2 = 0;
            uint64_t s3 = 0;
            for(long unsigned int i = 0; i < sizeof JUMP / sizeof *JUMP; i++)
                for(int b = 0; b < 64; b++) {
                    if (JUMP[i] & UINT64_C(1) << b) {
                        s0 ^= s[0];
                        s1 ^= s[1];
                        s2 ^= s[2];
                        s3 ^= s[3];
                    }
                    next();	
                }
                
            s[0] = s0;
            s[1] = s1;
            s[2] = s2;
            s[3] = s3;
        }


        /* This is the long-jump function for the generator. It is equivalent to
           2^192 calls to next(); it can be used to generate 2^64 starting points,
           from each of which jump() will generate 2^64 non-overlapping
           subsequences for parallel distributed computations. */

        void long_jump(void) {
            static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };

            uint64_t s0 = 0;
            uint64_t s1 = 0;
            uint64_t s2 = 0;
            uint64_t s3 = 0;
            for(long unsigned int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++)
                for(int b = 0; b < 64; b++) {
                    if (LONG_JUMP[i] & UINT64_C(1) << b) {
                        s0 ^= s[0];
                        s1 ^= s[1];
                        s2 ^= s[2];
                        s3 ^= s[3];
                    }
                    next();	
                }
                
            s[0] = s0;
            s[1] = s1;
            s[2] = s2;
            s[3] = s3;
        }
    };
};

