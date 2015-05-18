#include <util.h> // for discrete map
#include <random.h>

// GEN - INFINITE ROD - base estimators for a variety of step-length distributions and phase functions for multiple scattering in
//                    an infinite rod - positive x is right - source is at x = 0

class GENInfiniteRod
{
public: 

    virtual double samplestep() = 0; // step-length distribution sampling
    virtual double sampledir(const double prevDir ) = 0; // direction sampling
    virtual double init(double & pos, double & dir) = 0; // init particle position, weight, dir
    virtual void printDescriptor() = 0; // print a string describing the medium's variety

    // estimate albedo, internal distribution and moments for a finite rod scattering problem with vacuum boundary conditions
    void infiniteRodEstimatorAnalog(  const double c,           // single-scattering albedo
                                      const double maxx,        // maximum depth x in rod to estimate internal distribution
                                      const double dx,          // depth resolution for internal distribution
                                      const size_t numsamples, 
                                      const size_t numCollisionOrders,
                                      const size_t numMoments
                                  )
    {
        const size_t num_x_bins = floor( 2.0 * maxx / dx ) + 1.0;

        // radiance / vector flux
        size_t * collisionDensityL = new size_t [num_x_bins];
        size_t * collisionDensityR = new size_t [num_x_bins];

        for( size_t i = 0; i < num_x_bins; ++i )
        {
            collisionDensityL[i] = 0;
            collisionDensityR[i] = 0;
        }

        // collision density of arbitrary order
        const size_t num_collision_bins = num_x_bins * numCollisionOrders;
        size_t * nthcollisionDensitiesL = new size_t [num_collision_bins];
        size_t * nthcollisionDensitiesR = new size_t [num_collision_bins];
        for( size_t i = 0; i < num_collision_bins; ++i )
        {
            nthcollisionDensitiesL[i] = 0;
            nthcollisionDensitiesR[i] = 0;
        }

        // moments of fluence / density / scalar flux
        double * all_orders_moments = new double [numMoments];
        for( size_t i = 0; i < numMoments; ++i )
        {
            all_orders_moments[i] = 0;
        }
        double * per_order_moments = new double [numMoments * numCollisionOrders];
        for( size_t i = 0; i < numMoments * numCollisionOrders; ++i )
        {
            per_order_moments[i] = 0;
        }

        // perform samples
        for( size_t i = 0; i < numsamples; ++i )
        {
            size_t collisionOrder = 0;
            double pos( 0.0 );
            double dir( 1.0 );

            init( pos, dir );

            double step = samplestep();

            pos += dir * step;

            while( true )
            {
                collisionOrder++; // we collided - deal with absorption later

                const size_t x_index = discreteMap( -maxx, maxx, dx, pos );

                if( dir > 0 )
                {
                    collisionDensityR[x_index]++;
                }
                else
                {
                    collisionDensityL[x_index]++;
                }

                for( size_t m = 0; m < numMoments; ++m )
                {
                    all_orders_moments[m] += pow( pos, double( m ) );
                }

                if( collisionOrder < numCollisionOrders )
                {   
                    for( size_t m = 0; m < numMoments; ++m )
                    {
                        per_order_moments[m * numCollisionOrders + collisionOrder - 1] += pow( pos, double( m ) );
                    }
                }

                if( collisionOrder < numCollisionOrders )
                {   
                    const size_t collision_i = x_index + num_x_bins * ( collisionOrder - 1 );
                    
                    if( dir > 0 )
                    {
                        nthcollisionDensitiesR[collision_i]++;
                    }
                    else
                    {
                        nthcollisionDensitiesL[collision_i]++;
                    }
                }

                if( RandomReal() > c )
                {
                    // absorption - die
                    break;
                }

                dir = sampledir( dir );
                step = samplestep();
                pos += dir * step;
            }
        }

        printDescriptor();
        std::cout << "Analogmu_t c: " << c << 
                     " maxx: " << maxx << 
                     " dx: " << dx << 
                     " numsamples: " << numsamples << 
                     " numCollisionOrders: " << numCollisionOrders << 
                     " numMoments: " << numMoments << 
                     std::endl;
                     
        std::cout << "Collision density L:\n";
        for( size_t i = 0; i < num_x_bins; ++i )
        {
            std::cout << double( collisionDensityL[i] ) / ( double( numsamples ) * dx ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Collision density R:\n";
        for( size_t i = 0; i < num_x_bins; ++i )
        {
            std::cout << double( collisionDensityR[i] ) / ( double( numsamples ) * dx ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Density Moments:\n";
        for( size_t m = 0; m < numMoments; m++ )
        {
            std::cout << all_orders_moments[m] / ( double( numsamples ) ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Per-collision Density Moments:\n";
        for( size_t order = 0; order < numCollisionOrders; order++ )
        {
            for( size_t m = 0; m < numMoments; ++m )
            {
                std::cout << double( per_order_moments[order + m * numCollisionOrders ] ) / ( double( numsamples ) ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Per-collision Densities L:\n";
        for( size_t ci = 0; ci < numCollisionOrders; ++ci )
        {
            for( size_t xi = 0; xi < num_x_bins; ++xi )
            {
                const size_t collision_i = xi + num_x_bins * ci;
                std::cout << double( nthcollisionDensitiesL[collision_i] ) / ( double( numsamples ) * dx ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Per-collision Densities L:\n";
        for( size_t ci = 0; ci < numCollisionOrders; ++ci )
        {
            for( size_t xi = 0; xi < num_x_bins; ++xi )
            {
                const size_t collision_i = xi + num_x_bins * ci;
                std::cout << double( nthcollisionDensitiesR[collision_i] ) / ( double( numsamples ) * dx ) << " ";
            }
            std::cout << std::endl;
        }

        delete [] collisionDensityR;
        delete [] collisionDensityL;
        delete [] nthcollisionDensitiesL;
        delete [] nthcollisionDensitiesR;
    }
};
