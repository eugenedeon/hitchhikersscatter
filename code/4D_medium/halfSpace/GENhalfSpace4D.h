#include <util.h> // for discrete map
#include <random.h>

// GEN - 4D - HALF-SPACE - base estimator for a variety of step-length distributions and phase functions

class GENHalfspace
{
public: 

    virtual double samplestep() = 0; // step-length distribution sampling
    virtual Vector4 sampledir(const Vector4& prevDir ) = 0; // direction sampling
    virtual void printDescriptor() = 0; // print a string describing the medium's variety
    virtual Vector4 initPos() = 0;
    virtual Vector4 initDir() = 0;

    // estimate albedo, internal distribution and moments for a half-space scattering problem with vacuum boundary conditions
    void HalfSpaceEstimatorAnalog(  const double c,               // single-scattering albedo
                                    const double du,              // resolution of exitant angular distribution bins (u = cos theta )
                                    const double dz,              // resolution of depth-based quantities
                                    const double maxz,            // maximum depth recorded
                                    const size_t numsamples, 
                                    const size_t numCollisionOrders,
                                    const size_t numMoments
                                  )
    {
        const size_t numubins = size_t( floor( 1.0 / du ) ) + 1;
        const size_t numubins_internal = size_t( floor( 2.0 / du ) ) + 1;

        size_t countExitingParticles = 0; // used to estimate total albedo

        size_t * exitantDistribution = new size_t[numubins];
        size_t * exitantDistributionSingle = new size_t[numubins];
        size_t * exitantDistributionDouble = new size_t[numubins];
        for( size_t i = 0; i < numubins; ++i )
        {
            exitantDistribution[i] = 0;
            exitantDistributionSingle[i] = 0; // single-scattering
            exitantDistributionDouble[i] = 0; // double-scattering
        }

        size_t * perOrderExitingParticles = new size_t [numCollisionOrders];
        for( size_t i = 0; i < numCollisionOrders; ++i )
        {
            perOrderExitingParticles[i] = 0; // used to estimate per-collision order albedo 
        }      

        // internal density
        const size_t num_z_bins = size_t( floor( maxz / dz ) ) + 1;
        size_t * density = new size_t [num_z_bins];

        for( size_t i = 0; i <= num_z_bins; ++i )
        {
            density[i] = 0;
        }

        // internal angular collision density
        const size_t num_internal_angular_density_entries = num_z_bins * numubins_internal;
        size_t * density_angular = new size_t [num_internal_angular_density_entries];

        for( size_t i = 0; i <= num_internal_angular_density_entries; ++i )
        {
            density_angular[i] = 0;
        }

        // collision density of arbitrary order
        const size_t num_collision_bins = num_z_bins * numCollisionOrders;
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
            per_order_moments[i] = 0.0;
        }

        // angular moments
        double * angular_moments = new double[numMoments];
        for( size_t i = 0; i < numMoments; ++i )
        {
            angular_moments[i] = 0.0;
        }

        // perform samples
        for( size_t i = 0; i < numsamples; ++i )
        {
            size_t collisionOrder = 0;

            Vector4 pos( initPos() );
            Vector4 dir( initDir() );

            double step = samplestep();

            pos += dir * step;

            while( true )
            {    
                const double u = -dir.y;            
                if( pos.y <= 0.0 )
                {
                    //escape
                    countExitingParticles++;
                    if( collisionOrder < numCollisionOrders )
                    {   
                        perOrderExitingParticles[collisionOrder]++;
                    }
                    
                    const size_t u_index = discreteMap( 0.0, 1.0, du, u );
                    exitantDistribution[u_index]++;
                    if( 1 == collisionOrder ) exitantDistributionSingle[u_index]++;
                    if( 2 == collisionOrder ) exitantDistributionDouble[u_index]++;

                    for( size_t m = 0; m < numMoments; ++m )
                    {
                        angular_moments[m] += pow( u, double( m ) );
                    }

                    break;
                }

                collisionOrder++; // we collided - deal with absorption later

                const size_t z_index = discreteMap( 0.0, maxz, dz, pos.y );

                density[z_index]++;

                const size_t u_index_internal = discreteMap( -1.0, 1.0, du, u );
                assert(z_index * numubins_internal + u_index_internal < num_internal_angular_density_entries);
                density_angular[z_index * numubins_internal + u_index_internal]++;

                for( size_t m = 0; m < numMoments; ++m )
                {
                    all_orders_moments[m] += pow( pos.y, double( m ) );
                }

                if( collisionOrder < numCollisionOrders )
                {   
                    for( size_t m = 0; m < numMoments; ++m )
                    {
                        per_order_moments[m * numCollisionOrders + collisionOrder] += pow( pos.y, double( m ) );
                    }
                }

                if( collisionOrder < numCollisionOrders )
                {   
                    const size_t collision_i = z_index + num_z_bins * ( collisionOrder - 1 );
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
            }
        }

        printDescriptor();
        std::cout << "Analogmu_t c: " << c << 
                     " du: " << du << 
                     " maxz: " << maxz << 
                     " dz: " << dz << 
                     " numsamples: " << numsamples << 
                     " numCollisionOrders: " << numCollisionOrders << 
                     " numMoments: " << numMoments <<
                     std::endl;

        std::cout << "total exitance: " << double( countExitingParticles ) / double( numsamples ) << std::endl;
        std::cout << "per-collision order exitances: " << std::endl;
        for( size_t i = 0; i < numCollisionOrders; ++i )
        {
            std::cout << double( perOrderExitingParticles[i] ) / ( double( numsamples ) ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Exitant distribution: " << std::endl;
        for( size_t i = 0; i < numubins; ++i )
        {
            std::cout << double( exitantDistribution[i] ) / ( double( numsamples ) * du ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Exitant distribution single: " << std::endl;
        for( size_t i = 0; i < numubins; ++i )
        {
            std::cout << double( exitantDistributionSingle[i] ) / ( double( numsamples ) * du ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Exitant distribution double: " << std::endl;
        for( size_t i = 0; i < numubins; ++i )
        {
            std::cout << double( exitantDistributionDouble[i] ) / ( double( numsamples ) * du ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Internal density:\n";
        for( size_t i = 0; i < num_z_bins; ++i )
        {
            std::cout << double( density[i] ) / ( double( numsamples ) * dz ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Density Moments:\n";
        for( size_t m = 0; m < numMoments; m++ )
        {
            std::cout << all_orders_moments[m] / ( double( numsamples ) ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Angular Moments:\n";
        for( size_t m = 0; m < numMoments; m++ )
        {
            std::cout << angular_moments[m] / ( double( numsamples ) ) << " ";
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
            for( size_t ri = 0; ri < num_z_bins; ++ri )
            {
                const size_t collision_i = ri + num_z_bins * ( ci - 1 );
                std::cout << double( collisionDensities[collision_i] ) / ( double( numsamples ) * dz ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Internal Angular Collision Densities:\n";
        for( size_t zi = 0; zi < num_z_bins; ++zi )
        {
            for( size_t ui = 0; ui < numubins_internal; ++ui )
            {
                const size_t collision_angular_i = zi * numubins_internal + ui;
                std::cout << double( density_angular[collision_angular_i] ) / ( double( numsamples ) * 2.0 * du * dz ) << " ";
            }
            std::cout << std::endl;
        }

        delete [] density;
        delete [] density_angular;
    }
};
