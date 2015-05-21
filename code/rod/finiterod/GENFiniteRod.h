#include <util.h> // for discrete map
#include <random.h>

// GEN - FINITE ROD - base estimators for a variety of step-length distributions and phase functions for multiple scattering in
//                    a finite rod - positive x is right - boundaries are at x = 0 and x = a

class GENFiniteRod
{
public: 

    virtual double samplestep() = 0; // step-length distribution sampling
    virtual double sampledir(const double prevDir ) = 0; // direction sampling
    virtual double init(double & pos, double & dir) = 0; // init particle position, weight, dir
    virtual void printDescriptor() = 0; // print a string describing the medium's variety

    // estimate albedo, internal distribution and moments for a finite rod scattering problem with vacuum boundary conditions
    void RodAlbedoEstimatorAnalog(  const double c,           // single-scattering albedo
                                    const double rod_length,  // rod length
                                    const double dx,          // depth resolution for internal distribution
                                    const size_t numsamples, 
                                    const size_t numCollisionOrders,
                                    const size_t numMoments
                                  )
    {
        size_t countExitingParticlesR = 0; // used to estimate total reflectance albedo
        size_t countExitingParticlesT = 0; // used to estimate total transmitted albedo

        size_t * perOrderExitingParticlesR = new size_t [numCollisionOrders];
        size_t * perOrderExitingParticlesT = new size_t [numCollisionOrders];
        for( size_t i = 0; i < numCollisionOrders; ++i )
        {
            perOrderExitingParticlesR[i] = 0; // used to estimate per-collision order reflectance albedo 
            perOrderExitingParticlesT[i] = 0; // used to estimate per-collision order transmittance albedo 
        }

        const size_t num_x_bins = floor( rod_length / dx ) + 1.0;

        // left/right moving collision densities (direction of particles BEFORE collision)
        size_t * collisionDensityL = new size_t [num_x_bins];
        size_t * collisionDensityR = new size_t [num_x_bins];

        for( size_t i = 0; i < num_x_bins; ++i )
        {
            collisionDensityL[i] = 0;
            collisionDensityR[i] = 0;
        }

        // collision density of arbitrary order
        const size_t num_collision_bins = num_x_bins * numCollisionOrders;
        size_t * collisionDensities = new size_t [num_collision_bins];
        for( size_t i = 0; i < num_collision_bins; ++i )
        {
            collisionDensities[i] = 0;
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
                if ( pos >= rod_length )
                {
                    //escape
                    countExitingParticlesT++;
                    if( collisionOrder < numCollisionOrders )
                    {   
                        perOrderExitingParticlesT[collisionOrder]++;
                    }
                    break;
                }
                collisionOrder++; // we collided - deal with absorption later

                const size_t x_index = discreteMap( 0.0, rod_length, dx, pos );

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
                        per_order_moments[m * numCollisionOrders + collisionOrder] += pow( pos, double( m ) );
                    }
                }

                if( collisionOrder < numCollisionOrders )
                {   
                    const size_t collision_i = x_index + num_x_bins * ( collisionOrder - 1 );
                    collisionDensities[collision_i]++;
                }

                if( RandomReal() > c )
                {
                    // absorption - die
                    break;
                }

                dir = sampledir( dir );
                step = samplestep();
                pos += dir * step;

                if( pos <= 0.0 )
                {
                    //escape
                    countExitingParticlesR++;
                    if( collisionOrder < numCollisionOrders )
                    {   
                        perOrderExitingParticlesR[collisionOrder]++;
                    }
                    break;
                }
            }
        }

        printDescriptor();
        std::cout << "Analogmu_t c: " << c <<
                     " rod_length: " << rod_length << 
                     " dx: " << dx << 
                     " numsamples: " << numsamples << 
                     " numCollisionOrders: " << numCollisionOrders << 
                     " numMoments: " << numMoments << 
                     std::endl;

        std::cout << "total albedo R: " << double( countExitingParticlesR ) / double( numsamples ) << std::endl;
        std::cout << "total albedo T: " << double( countExitingParticlesT ) / double( numsamples ) << std::endl;
        
        std::cout << "per-collision order albedos R: " << std::endl;
        for( size_t i = 0; i < numCollisionOrders; ++i )
        {
            std::cout << double( perOrderExitingParticlesR[i] ) / ( double( numsamples ) ) << " ";
        }
        std::cout << std::endl;
        
        std::cout << "per-collision order albedos T: " << std::endl;
        for( size_t i = 0; i < numCollisionOrders; ++i )
        {
            std::cout << double( perOrderExitingParticlesT[i] ) / ( double( numsamples ) ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Left-moving collision density:\n";
        for( size_t i = 0; i < num_x_bins; ++i )
        {
            std::cout << double( collisionDensityL[i] ) / ( double( numsamples ) * dx ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Right-moving collision density:\n";
        for( size_t i = 0; i < num_x_bins; ++i )
        {
            std::cout << double( collisionDensityR[i] ) / ( double( numsamples ) * dx ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Collision Density Moments:\n";
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

        std::cout << "Per-collision Densities:\n";
        for( size_t ci = 0; ci < numCollisionOrders; ++ci )
        {
            for( size_t ri = 0; ri < num_x_bins; ++ri )
            {
                const size_t collision_i = ri + num_x_bins * ( ci - 1 );
                std::cout << double( collisionDensities[collision_i] ) / ( double( numsamples ) * dx ) << " ";
            }
            std::cout << std::endl;
        }

        delete [] collisionDensityR;
        delete [] collisionDensityL;
    }
};

