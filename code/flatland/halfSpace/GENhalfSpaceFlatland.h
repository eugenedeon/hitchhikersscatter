#include <util.h> // for discrete map
#include <random.h>

// GEN - FLATLAND - HALF-SPACE - base estimator for a variety of step-length distributions and phase functions

class GENHalfspace
{
public: 

    virtual double samplestep() = 0; // step-length distribution sampling
    virtual Vector2 sampledir(const Vector2& prevDir ) = 0; // direction sampling
    virtual double mfp() = 0;
    virtual void printDescriptor() = 0; // print a string describing the medium's variety
    virtual Vector2 initPos() = 0;
    virtual Vector2 initDir() = 0;

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

        double track_length = 0.0;
        double expected_track_length = 0.0;
        double expected_track_length2 = 0.0;

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

            Vector2 pos( initPos() );
            Vector2 dir( initDir() );

            double step = samplestep();
            track_length += step;

            // expected track length upon entry is simple the uncorrelated mfp
            expected_track_length += mfp();
            expected_track_length2 += mfp();

            pos += dir * step;

            while( true )
            {    
                const double u = -dir.y;            
                if( pos.y <= 0.0 )
                {
                    //escape

                    // how much did we over escape by?
                    track_length -= pos.y / dir.y;

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
                track_length += step; // too much if we escape here - subtract later

                expected_track_length += ( step <= pos.y ? step : (M_PI*step - step*acos(pos.y/step) + pos.y*
      log((step + sqrt(-pos.y*pos.y + step*step))/pos.y))/M_PI );
                expected_track_length2 += ( step <= pos.y ? step : (M_PI*step - step*acos(pos.y/step) + pos.y*
      log((step + sqrt(-pos.y*pos.y + step*step))/pos.y))/M_PI ) - 0.5 * step + 0.5 * mfp();

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

        std::cout << "density moment 0: " << all_orders_moments[0] / ( double( numsamples ) ) << std::endl;
        std::cout << "Track length: " << track_length / double( numsamples ) << std::endl;
        std::cout << "Expected track length: " << expected_track_length / double( numsamples ) << std::endl;
        std::cout << "Expected track length 2: " << expected_track_length2 / double( numsamples ) << std::endl;

        delete [] density;
        delete [] density_angular;
    }

    // estimate total exitance flux at radii r from surface position of normally-incident searchlight illumination
    // also output laterally-integrtaed angular exitance (in the BRDF sense)
    void searchLightNormalIncidenceAnalog(   const double c, // single-scattering albedo
                                             const double maxr, // maximum radius to record
                                             const double dr, // uniformly spaced radial bins of with dr
                                             const double du, // uniformly spaced direction cosine bins
                                             const size_t numsamples
                                          )
    {
        size_t total_exitance_count = 0;

        const size_t num_angular_bins = floor( 1.0 / du ) + 1;
        const size_t num_radial_bins = floor( maxr / dr ) + 1;
        
        size_t * angular_exitance = new size_t [num_angular_bins];
        size_t * radial_exitance = new size_t [num_radial_bins];
        
        for( size_t i = 0; i < num_angular_bins; ++i )
        {
            angular_exitance[i] = 0;
        }

        for( size_t i = 0; i < num_radial_bins; ++i )
        {
            radial_exitance[i] = 0;
        }

        for( size_t i = 0; i < numsamples; ++i )
        {
            Vector2 pos( 0.0, 0.0 );
            Vector2 dir( 0.0, 1.0 ); // 2D searchlight: normal incidence at the origin

            double step = samplestep();

            pos += dir * step;

            while( RandomReal() < c ) // absorption = death
            {
                dir = sampledir(dir);
                step = samplestep();
                pos += dir * step;

                if( pos.y <= 0.0 )
                {
                    // escape and score
                    size_t angular_index = discreteMap( 0.0, 1.0, du, -dir.y );
                    angular_exitance[ angular_index ]++;

                    // move back to the surface
                    const double tbackup = -pos.y / dir.y;
                    pos += dir * tbackup;

                    const double r = Norm(pos);
                    size_t radial_index = discreteMap( 0.0, maxr, dr, r );
                    radial_exitance[ radial_index ]++;

                    total_exitance_count++;
                    break;
                }
            }
        }

        printDescriptor();
        std::cout << "Analogmu_t c: " << c << 
                     " maxr: " << maxr << 
                     " dr: " << dr << 
                     " du: " << du << 
                     " numsamples: " << numsamples << 
                     std::endl;

        std::cout << "global albedo: " << double( total_exitance_count ) / double( numsamples ) << std::endl;

        std::cout << "angular exitance:\n";
        for( size_t i = 0; i < num_angular_bins; ++i )
        {
            std::cout << ( angular_exitance[i] ) / double( numsamples * du ) << " ";
        }

        std::cout << std::endl;

        std::cout << "radial exitant flux:\n";
        for( size_t i = 0; i < num_radial_bins; ++i )
        {
            std::cout << ( radial_exitance[i] ) / double( numsamples * dr ) << " ";
        }

        std::cout << std::endl;

        delete [] angular_exitance;
        delete [] radial_exitance;

    };
};
