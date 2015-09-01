#include <util.h> // for discrete map
#include <random.h>

// GEN - FLATLAND - INFINITE MEDIUM - base estimator for a variety of step-length distributions and phase functions

class GENInfiniteMedium
{
public: 

    virtual double samplestep() = 0; // step-length distribution sampling
    virtual Vector2 sampledir(const Vector2& prevDir ) = 0; // direction sampling
    virtual void printDescriptor() = 0; // print a string describing the medium's variety

    // estimate fluence at various radii from an isotropic point source
    void isotropicPointSourceAnalogMut(  const double c,
                                         const double mu_t, 
                                         const double maxr, 
                                         const double dr, // uniformly spaced radial bins
                                         const double du, // uniformly spaced direction cosine bins
                                         const size_t numsamples, 
                                         const size_t numCollisionOrders,
                                         const size_t numMoments
                                      )
    {
        const double mu_a = mu_t * ( 1.0 - c );
        const double mu_s = mu_t * c;

        // collision density
        const size_t num_r_bins = floor( maxr / dr ) + 1.0;
        size_t * collisionDensity = new size_t [num_r_bins];
        for( size_t i = 0; i < num_r_bins; ++i )
        {
            collisionDensity[i] = 0;
        }

        // angular collision density
        const size_t num_u_bins = floor( 2.0 / du ) + 1.0;
        const size_t num_ru_bins = num_r_bins * num_u_bins;
        size_t * angularCollisionDensity = new size_t [num_ru_bins];
        for( size_t i = 0; i < num_ru_bins; ++i )
        {
            angularCollisionDensity[i] = 0;
        }

        // collision density of arbitrary order
        const size_t numCollisionbins = num_r_bins * numCollisionOrders;
        size_t * collisionDensities = new size_t [numCollisionbins];
        for( size_t i = 0; i < numCollisionbins; ++i )
        {
            collisionDensities[i] = 0;
        }

        // moments
        double * globalMoments = new double [numMoments];
        for( size_t i = 0; i < numMoments; ++i )
        {
            globalMoments[i] = 0;
        }

        const size_t numOrderMoments = numMoments * numCollisionOrders;
        double * orderMoments = new double [numOrderMoments];
        for( size_t i = 0; i < numOrderMoments; ++i )
        {
            orderMoments[i] = 0;
        }

        for( size_t i = 0; i < numsamples; ++i )
        {
            size_t collisionOrder = 0;

            Vector2 pos( 0.0, 0.0 );
            Vector2 dir( 1.0, 0.0 );

            double step = samplestep();

            pos += dir * step;

            while( RandomReal() < c )
            {
                const double r = Norm( pos );
                collisionOrder++; // we collided and didn't absorb

                for( size_t m = 0; m < numMoments; ++m )
                {
                    globalMoments[m] += pow( r, double( m ) );
                }

                if( collisionOrder < numCollisionOrders )
                {   
                    for( size_t m = 0; m < numMoments; ++m )
                    {
                        assert( m * numCollisionOrders + collisionOrder < numOrderMoments );
                        orderMoments[m * numCollisionOrders + collisionOrder] += pow( r, double( m ) );
                    }
                }

                const size_t ri = std::min( num_r_bins - 1, (long unsigned int)( floor( r / dr ) ) );
                assert( ri < num_r_bins );
                collisionDensity[ri]++;

                const float u = Dot( Normalize(pos), dir );
                const size_t ui = discreteMap( -1.0, 1.0, du, u );
                const size_t rui = num_u_bins * ri + ui;
                assert( rui < num_ru_bins );
                angularCollisionDensity[rui]++;

                if( collisionOrder < numCollisionOrders )
                {   
                    const size_t collision_i = ri + num_r_bins * ( collisionOrder - 1 );
                    assert( collision_i < numCollisionbins );
                    collisionDensities[collision_i]++;
                }

                dir = sampledir( dir );
                step = samplestep();
                pos += dir * step;
            }

            collisionOrder++; // we collided and were absorbed

            const double r = Norm( pos );
            for( size_t m = 0; m < numMoments; ++m )
            {
                globalMoments[m] += pow( r, double( m ) );
            }

            if( collisionOrder < numCollisionOrders )
            {   
                for( size_t m = 0; m < numMoments; ++m )
                {
                    assert( m * numCollisionOrders + collisionOrder < numOrderMoments );
                    orderMoments[m * numCollisionOrders + collisionOrder] += pow( r, double( m ) );
                }
            }
            
            const size_t ri = std::min( num_r_bins - 1, (long unsigned int)( floor( r / dr ) ) );
            assert( ri < num_r_bins );
            collisionDensity[ri]++;

            const float u = Dot( Normalize(pos), dir );
            const size_t ui = discreteMap( -1.0, 1.0, du, u );
            const size_t rui = num_u_bins * ri + ui;
            assert( rui < num_ru_bins );
            angularCollisionDensity[rui]++;

            if( collisionOrder < numCollisionOrders )
            {   
                const size_t collision_i = ri + num_r_bins * ( collisionOrder - 1 );
                assert( collision_i < numCollisionbins );
                collisionDensities[collision_i]++;
            }
        }

        printDescriptor();
        std::cout << "Analogmu_t c: " << c << 
                     " mu_t: " << mu_t << 
                     " maxr: " << maxr << 
                     " dr: " << dr << 
                     " du: " << du << 
                     " numsamples: " << numsamples << 
                     " numCollisionOrders: " << numCollisionOrders << 
                     " numMoments: " << numMoments << 
                     std::endl;

        std::cout << "CollisionDensity:\n";
        for( size_t i = 0; i < num_r_bins; ++i )
        {
            std::cout << double( collisionDensity[i] ) / ( double( numsamples ) * dr ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Density Moments:\n";
        for( size_t m = 0; m < numMoments; m++ )
        {
            std::cout << globalMoments[m] / double( numsamples ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Per-collision Density Moments:\n";
        for( size_t order = 0; order < numCollisionOrders; order++ )
        {
            for( size_t m = 0; m < numMoments; ++m )
            {
                std::cout << double( orderMoments[order + m * numCollisionOrders ] ) / double( numsamples ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Per-collision Densities:\n";
        for( size_t ci = 0; ci < numCollisionOrders; ++ci )
        {
            for( size_t ri = 0; ri < num_r_bins; ++ri )
            {
                const size_t collision_i = ri + num_r_bins * ( ci - 1 );
                std::cout << double( collisionDensities[collision_i] ) / ( double( numsamples ) * dr ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Angular Collision Densities:\n";
        for( size_t ri = 0; ri < num_r_bins; ++ri )
        {
            for( size_t ui = 0; ui < num_u_bins; ++ui )
            {
                const size_t rui = ri * num_u_bins + ui;
                std::cout << double( angularCollisionDensity[rui] ) / ( double( numsamples ) * dr * du ) << " ";
            }
            std::cout << std::endl;
        }

        delete [] collisionDensity;
    }
};
