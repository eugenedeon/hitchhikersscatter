#include <util.h> // for discrete map
#include <random.h>

// GEN - 3D - INFINITE MEDIUM GRT - base estimator for a variety of step-length distributions and phase functions

class GENInfiniteMedium
{
public: 

    virtual double sample_correlated_step() = 0; // correlated step-length distribution sampling
    virtual double sample_uncorrelated_step() = 0; // uncorrelated step-length distribution sampling

    virtual Vector3 sampledir(const Vector3& prevDir ) = 0; // direction sampling

    // step-length-dependent cross sections and their integrals for collision and track-length estimators
    virtual double sigma_tc( const double s ) = 0;
    virtual double sigma_tu( const double s ) = 0;
    virtual double correlated_collision_integral( const double s1, const double s2 ) = 0;
    virtual double uncorrelated_collision_integral( const double s1, const double s2 ) = 0;
    
    virtual void printDescriptor() = 0; // print a string describing the medium's variety

    // estimate quantities at various radii from an isotropic point source - Analog walk, collision estimators of flux and reaction rate
    void isotropicPointSourceAnalogCollisionEstimator(  const double c, // single-scattering albedo
                                                        const double maxr, // maximum radius to record
                                                        const double dr, // uniformly spaced radial bins of with dr
                                                        const double maxt,
                                                        const double dt,
                                                        const double du, // uniformly spaced direction cosine bins
                                                        const size_t numsamples, 
                                                        const size_t numCollisionOrders,
                                                        const size_t numMoments
                                                     )
    {
        // scalar quantities: 
        //   collisionDensity - collision rate density for entering collisions at distance r from point source
        //   fluence
        const size_t num_r_bins = floor( maxr / dr ) + 1.0;
        size_t * collisionDensity = new size_t [num_r_bins];
        double * fluence = new double [num_r_bins];
        for( size_t i = 0; i < num_r_bins; ++i )
        {
            collisionDensity[i] = 0;
            fluence[i] = 0.0;
        }

        // time-resolved collision rate density and fluence tallies
        const size_t num_t_bins = floor( maxt / dt ) + 1.0;
        const size_t num_tr_bins = num_r_bins * num_t_bins;
        size_t * collisionDensityTR = new size_t [num_tr_bins];
        double * fluenceTR = new double [num_tr_bins];
        for( size_t i = 0; i < num_tr_bins; ++i )
        {
            collisionDensityTR[i] = 0;
            fluenceTR[i] = 0.0;
        }

        // angular quantities:
        //  angularCollisionDensity - collision rate density for ENTERING collisions at distance r from point source 
        //                            along direction with cosine u to the position vector
        //  angularSourceDensity    - rate density of particles LEAVING a collision at distance r from point source
        //                            into direction with cosine u to the position vector
        //  angularFlux             - azimuthally-integrated radiance
        const size_t num_u_bins = floor( 2.0 / du ) + 1.0;
        const size_t num_ru_bins = num_r_bins * num_u_bins;
        size_t * angularCollisionDensity = new size_t [num_ru_bins];
        size_t * angularSourceDensity = new size_t [num_ru_bins];
        double * angularFlux = new double [num_ru_bins];
        for( size_t i = 0; i < num_ru_bins; ++i )
        {
            angularCollisionDensity[i] = 0;
            angularSourceDensity[i] = 0;
            angularFlux[i] = 0.0;
        }

        // scalar quantities by collision order
        const size_t numCollisionbins = num_r_bins * numCollisionOrders;
        size_t * collisionDensities = new size_t [numCollisionbins];
        double * fluences = new double [numCollisionbins];
        for( size_t i = 0; i < numCollisionbins; ++i )
        {
            collisionDensities[i] = 0;
            fluences[i] = 0.0;
        }

        // spatial moments
        double * spatialCollisionRateMoments = new double [numMoments];
        double * spatialFluenceMoments = new double [numMoments];
        for( size_t i = 0; i < numMoments; ++i )
        {
            spatialCollisionRateMoments[i] = 0;
            spatialFluenceMoments[i] = 0;
        }

        // spatial moments by collision order
        const size_t numOrderMoments = numMoments * numCollisionOrders;
        double * spatialCollisionRateMomentsByCollision = new double [numOrderMoments];
        double * spatialFluenceMomentsByCollision = new double [numOrderMoments];
        for( size_t i = 0; i < numOrderMoments; ++i )
        {
            spatialCollisionRateMomentsByCollision[i] = 0;
            spatialFluenceMomentsByCollision[i] = 0;
        }

        double track_length(0.0);

        for( size_t i = 0; i < numsamples; ++i )
        {
            size_t collisionOrder = 0;
            double t = 0; // track-length of the current particle only
            Vector3 pos( 0.0, 0.0, 0.0 );
            Vector3 dir( 1.0, 0.0, 0.0 );
            bool absorbed = false;

            do
            {
                double step = ( 0 == collisionOrder ) ? sample_uncorrelated_step() : sample_correlated_step();
                pos += dir * step;
                track_length += step;
                t += step;

                collisionOrder++; // we collided
                const double r = Norm( pos );
                const size_t ri = discreteMap( 0.0, maxr, dr, r );
                const size_t ti = discreteMap( 0.0, maxt, dt, t );

                // collision estimator for angular flux/fluence
                const double fluxCollisionEstimatorScore = ( 1 == collisionOrder ) ? ( 1.0 / sigma_tu( step ) ) : ( 1.0 / sigma_tc( step ) );

                // tally scalar densities
                collisionDensity[ri]++;
                collisionDensityTR[ ri + num_r_bins * ti ]++;
                fluence[ri] += fluxCollisionEstimatorScore;
                fluenceTR[ ri + num_r_bins * ti ] += fluxCollisionEstimatorScore;

                // tally angular densities
                const float u = Dot( Normalize(pos), dir );
                const size_t ui = discreteMap( -1.0, 1.0, du, u );
                const size_t rui = num_u_bins * ri + ui;
                assert( rui < num_ru_bins );
                angularCollisionDensity[rui]++;
                angularFlux[rui] += fluxCollisionEstimatorScore;

                // tally spatial/radial moments
                for( size_t m = 0; m < numMoments; ++m )
                {
                    spatialCollisionRateMoments[m]  += pow( r, double( m ) );
                    spatialFluenceMoments[m]        += pow( r, double( m ) ) * fluxCollisionEstimatorScore;
                }

                if( collisionOrder < numCollisionOrders )
                {   
                    for( size_t m = 0; m < numMoments; ++m )
                    {
                        const size_t orderMomentIndex = m * numCollisionOrders + collisionOrder;
                        assert( orderMomentIndex < numOrderMoments );
                        spatialCollisionRateMomentsByCollision[orderMomentIndex] += pow( r, double( m ) );
                        spatialFluenceMomentsByCollision[orderMomentIndex - 1]   += pow( r, double( m ) ) * fluxCollisionEstimatorScore;
                    }
                }

                if( collisionOrder < numCollisionOrders )
                {   
                    const int collision_i = ri + num_r_bins * ( collisionOrder );
                    const int fluence_i = ri + num_r_bins * ( collisionOrder - 1 );
                    assert( 0 <= collision_i && collision_i < numCollisionbins );
                    assert( 0 <= fluence_i && fluence_i < numCollisionbins );
                    collisionDensities[collision_i]++;
                    fluences[fluence_i] += fluxCollisionEstimatorScore;
                }

                absorbed = RandomReal() > c;

                if( !absorbed )
                {
                    dir = sampledir( dir );
                    // score angular source density
                    const float u = Dot( Normalize(pos), dir );
                    const size_t ui = discreteMap( -1.0, 1.0, du, u );
                    const size_t rui = num_u_bins * ri + ui;
                    assert( rui < num_ru_bins );
                    angularSourceDensity[rui]++;
                }
            }
            while( !absorbed );
        }

        printDescriptor();
        std::cout << "AnalogCollisionEstimator c: " << c << 
                     " maxr: " << maxr << 
                     " dr: " << dr << 
                     " du: " << du << 
                     " numsamples: " << numsamples << 
                     " numCollisionOrders: " << numCollisionOrders << 
                     " numMoments: " << numMoments << 
                     std::endl;

        std::cout << "Collision-rate density for entering collisions at distance r from point source:\n";
        for( size_t i = 0; i < num_r_bins; ++i )
        {
            std::cout << double( collisionDensity[i] ) / ( double( numsamples ) * dr ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Fluence / scalar flux at distance r from point source:\n";
        for( size_t i = 0; i < num_r_bins; ++i )
        {
            std::cout << double( fluence[i] ) / ( double( numsamples ) * dr ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Spatial radial moments of the Collision-Rate Density:\n";
        for( size_t m = 0; m < numMoments; m++ )
        {
            std::cout << spatialCollisionRateMoments[m] / double( numsamples ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Spatial radial moments of the fluence:\n";
        for( size_t m = 0; m < numMoments; m++ )
        {
            std::cout << spatialFluenceMoments[m] / double( numsamples ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Per-collision Collision Density Moments:\n";
        for( size_t order = 0; order < numCollisionOrders; order++ )
        {
            for( size_t m = 0; m < numMoments; ++m )
            {
                std::cout << double( spatialCollisionRateMomentsByCollision[order + m * numCollisionOrders ] ) / double( numsamples ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Per-collision Fluence Moments:\n";
        for( size_t order = 0; order < numCollisionOrders; order++ )
        {
            for( size_t m = 0; m < numMoments; ++m )
            {
                std::cout << double( spatialFluenceMomentsByCollision[order + m * numCollisionOrders ] ) / double( numsamples ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Per-collision Collision Rate Densities:\n";
        for( int ci = 0; ci < numCollisionOrders; ++ci )
        {
            for( int ri = 0; ri < num_r_bins; ++ri )
            {
                const int collision_i = ri + num_r_bins * ci;
                assert( 0 <= collision_i && collision_i < numCollisionbins );
                std::cout << double( collisionDensities[collision_i] ) / ( double( numsamples ) * dr ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Per-collision Fluences:\n";
        for( int ci = 0; ci < numCollisionOrders; ++ci )
        {
            for( int ri = 0; ri < num_r_bins; ++ri )
            {
                const int collision_i = ri + num_r_bins * ci;
                assert( 0 <= collision_i && collision_i < numCollisionbins );
                std::cout << double( fluences[collision_i] ) / ( double( numsamples ) * dr ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Angular Collision Rate Densities:\n";
        for( size_t ri = 0; ri < num_r_bins; ++ri )
        {
            for( size_t ui = 0; ui < num_u_bins; ++ui )
            {
                const size_t rui = ri * num_u_bins + ui;
                std::cout << double( angularCollisionDensity[rui] ) / ( double( numsamples ) * dr * du ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Angular Flux:\n";
        for( size_t ri = 0; ri < num_r_bins; ++ri )
        {
            for( size_t ui = 0; ui < num_u_bins; ++ui )
            {
                const size_t rui = ri * num_u_bins + ui;
                std::cout << double( angularFlux[rui] ) / ( double( numsamples ) * dr * du ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Angular Collision Source Density:\n";
        for( size_t ri = 0; ri < num_r_bins; ++ri )
        {
            for( size_t ui = 0; ui < num_u_bins; ++ui )
            {
                const size_t rui = ri * num_u_bins + ui;
                std::cout << double( angularSourceDensity[rui] ) / ( double( numsamples ) * dr * du ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Time-resolved collision rate densities:\n";
        for( int ti = 0; ti < num_t_bins; ++ti )
        {
            for( int ri = 0; ri < num_r_bins; ++ri )
            {
                const int collision_i = ri + num_r_bins * ti;
                std::cout << double( collisionDensityTR[collision_i] ) / ( double( numsamples ) * dr * dt ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Time-resolved fluence:\n";
        for( int ti = 0; ti < num_t_bins; ++ti )
        {
            for( int ri = 0; ri < num_r_bins; ++ri )
            {
                const int collision_i = ri + num_r_bins * ti;
                std::cout << double( fluenceTR[collision_i] ) / ( double( numsamples ) * dr * dt ) << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "Mean track length: " << track_length / double( numsamples ) << "\n";
    }

    // estimate quantities at various depths from an isotropic plane source - Analog walk, collision estimators of flux and reaction rate
    void isotropicPlaneSourceAnalogCollision(  const double c, // single-scattering albedo
                                         const double maxr, // maximum radius to record
                                         const double dr, // uniformly spaced radial bins of with dr
                                         const double du, // uniformly spaced direction cosine bins
                                         const size_t numsamples, 
                                         const size_t numCollisionOrders,
                                         const size_t numMoments
                                      )
    {
        // collision density
        const size_t num_r_bins = floor( maxr / dr ) + 1.0;
        size_t * collisionDensity = new size_t [num_r_bins];
        double * fluence = new double [num_r_bins];
        for( size_t i = 0; i < num_r_bins; ++i )
        {
            collisionDensity[i] = 0;
            fluence[i] = 0.0;
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

        double * globalphiMoments = new double [numMoments];
        for( size_t i = 0; i < numMoments; ++i )
        {
            globalphiMoments[i] = 0;
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

            Vector3 pos( 0.0, 0.0, 0.0 );
            Vector3 dir( 1.0, 0.0, 0.0 );
            double step = sample_uncorrelated_step();

            pos += dir * step;

            while( RandomReal() < c )
            {
                const double r = pos.x;
                collisionOrder++; // we collided and didn't absorb

                for( size_t m = 0; m < numMoments; ++m )
                {
                    globalMoments[m] += pow( r, double( m ) );
                    globalphiMoments[m] += pow( r, double( m ) ) * (1 == collisionOrder ? 1.0 / sigma_tu( step ) : 1.0 / sigma_tc( step ));
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
                fluence[ri] += 1 == collisionOrder ? 1.0 / sigma_tu( step ) : 1.0 / sigma_tc( step );

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
                step = sample_correlated_step();
                pos += dir * step;
            }

            collisionOrder++; // we collided and were absorbed

            const double r = pos.x;
            for( size_t m = 0; m < numMoments; ++m )
            {
                globalMoments[m] += pow( r, double( m ) );
                globalphiMoments[m] += pow( r, double( m ) ) * (1 == collisionOrder ? 1.0 / sigma_tu( step ) : 1.0 / sigma_tc( step ));
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
            fluence[ri] += 1 == collisionOrder ? 1.0 / sigma_tu( step ) : 1.0 / sigma_tc( step );

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

        std::cout << "Flux Moments:\n";
        for( size_t m = 0; m < numMoments; m++ )
        {
            std::cout << globalphiMoments[m] / double( numsamples ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Fluence:\n";
        for( size_t i = 0; i < num_r_bins; ++i )
        {
            std::cout << double( fluence[i] ) / ( double( numsamples ) * dr ) << " ";
        }
        std::cout << std::endl;

        delete [] collisionDensity;
        delete [] fluence;
    }

    void isotropicPlaneSourceAnalogTrackLength(  const double c, // single-scattering albedo
                                                 const double maxz, // maximum spatial extent to tally
                                                 const double dz, // uniformly spaced spatial bins of thickness dz
                                                 const double du, // uniformly spaced direction cosine bins
                                                 const size_t numsamples, 
                                                 const size_t numCollisionOrders,
                                                 const size_t numMoments
                                              )
    {
        // spatial tallies
        const size_t num_spatial_bins = discreteMap( -maxz, maxz, dz, maxz ) + 1;
        double * collisionDensity = new double [num_spatial_bins];
        double * fluence = new double [num_spatial_bins];
        for( size_t i = 0; i < num_spatial_bins; ++i )
        {
            collisionDensity[i] = 0.0;
            fluence[i] = 0.0;
        }

        for( size_t i = 0; i < numsamples; ++i )
        {
            size_t collisionOrder = 0;

            Vector3 pos( 0.0, 0.0, 0.0 );
            Vector3 prevpos = pos;
            Vector3 dir( isotropicDir() );
            double step = sample_uncorrelated_step();
            pos += dir * step;

            for( size_t j = 0; j < num_spatial_bins - 1; ++j )
            {
                const double ztallymin = -maxz + dz * j;
                const double ztallymax = ztallymin + dz;

                // complete crossing
                if( ( pos.z < ztallymin && prevpos.z > ztallymax ) )
                {
                    fluence[j] += dz / fabs( dir.z );
                    collisionDensity[j] += uncorrelated_collision_integral( ( prevpos.z - ztallymax ) / fabs( dir.z ), ( prevpos.z - ztallymax + dz ) / fabs( dir.z ) );
                }

                // complete crossing - other direction
                if( ( prevpos.z < ztallymin && pos.z > ztallymax ) )
                {
                    fluence[j] += dz / fabs( dir.z );
                    collisionDensity[j] += uncorrelated_collision_integral( ( ztallymin - prevpos.z ) / fabs( dir.z ), ( ztallymin - prevpos.z + dz ) / fabs( dir.z ) );
                }

                // track is entirely contained inside tally region
                if( pos.z >= ztallymin && pos.z <= ztallymax && prevpos.z >= ztallymin && prevpos.z <= ztallymax )
                {
                    fluence[j] += step;
                    collisionDensity[j] += uncorrelated_collision_integral( 0.0, step );
                }

                if( prevpos.z >= ztallymin && prevpos.z <= ztallymax && pos.z > ztallymax )
                {
                    fluence[j] += ( ztallymax - prevpos.z ) / fabs( dir.z );
                    collisionDensity[j] += uncorrelated_collision_integral( 0.0, ( ztallymax - prevpos.z ) / fabs( dir.z ) );
                }

                if( pos.z >= ztallymin && pos.z <= ztallymax && prevpos.z > ztallymax )
                {
                    fluence[j] += ( ztallymax - pos.z ) / fabs( dir.z );
                    collisionDensity[j] += uncorrelated_collision_integral( ( prevpos.z - ztallymax ) / fabs( dir.z ), step );
                }

                if( prevpos.z >= ztallymin && prevpos.z <= ztallymax && pos.z < ztallymin )
                {
                    fluence[j] += ( prevpos.z - ztallymin ) / fabs( dir.z );
                    collisionDensity[j] += uncorrelated_collision_integral( 0.0, ( prevpos.z - ztallymin ) / fabs( dir.z ) );
                }

                if( pos.z >= ztallymin && pos.z <= ztallymax && prevpos.z < ztallymin )
                {
                    fluence[j] += ( pos.z - ztallymin ) / fabs( dir.z );
                    collisionDensity[j] += uncorrelated_collision_integral( ( ztallymin - prevpos.z ) / fabs( dir.z ), step );
                }
            }

            while( RandomReal() < c )
            {
                dir = sampledir( dir );
                step = sample_correlated_step();
                prevpos = pos;
                pos += dir * step;

                for( size_t j = 0; j < num_spatial_bins - 1; ++j )
                {
                    const double ztallymin = -maxz + dz * j;
                    const double ztallymax = ztallymin + dz;

                    // complete crossing
                    if( ( pos.z < ztallymin && prevpos.z > ztallymax ) )
                    {
                        fluence[j] += dz / fabs( dir.z );
                        collisionDensity[j] += correlated_collision_integral( ( prevpos.z - ztallymax ) / fabs( dir.z ), ( prevpos.z - ztallymax + dz ) / fabs( dir.z ) );
                    }

                    // complete crossing - other direction
                    if( ( prevpos.z < ztallymin && pos.z > ztallymax ) )
                    {
                        fluence[j] += dz / fabs( dir.z );
                        collisionDensity[j] += correlated_collision_integral( ( ztallymin - prevpos.z ) / fabs( dir.z ), ( ztallymin - prevpos.z + dz ) / fabs( dir.z ) );
                    }

                    // track is entirely contained inside tally region
                    if( pos.z >= ztallymin && pos.z <= ztallymax && prevpos.z >= ztallymin && prevpos.z <= ztallymax )
                    {
                        fluence[j] += step;
                        collisionDensity[j] += correlated_collision_integral( 0.0, step );
                    }

                    if( prevpos.z >= ztallymin && prevpos.z <= ztallymax && pos.z > ztallymax )
                    {
                        fluence[j] += ( ztallymax - prevpos.z ) / fabs( dir.z );
                        collisionDensity[j] += correlated_collision_integral( 0.0, ( ztallymax - prevpos.z ) / fabs( dir.z ) );
                    }

                    if( pos.z >= ztallymin && pos.z <= ztallymax && prevpos.z > ztallymax )
                    {
                        fluence[j] += ( ztallymax - pos.z ) / fabs( dir.z );
                        collisionDensity[j] += correlated_collision_integral( ( prevpos.z - ztallymax ) / fabs( dir.z ), step );
                    }

                    if( prevpos.z >= ztallymin && prevpos.z <= ztallymax && pos.z < ztallymin )
                    {
                        fluence[j] += ( prevpos.z - ztallymin ) / fabs( dir.z );
                        collisionDensity[j] += correlated_collision_integral( 0.0, ( prevpos.z - ztallymin ) / fabs( dir.z ) );
                    }

                    if( pos.z >= ztallymin && pos.z <= ztallymax && prevpos.z < ztallymin )
                    {
                        fluence[j] += ( pos.z - ztallymin ) / fabs( dir.z );
                        collisionDensity[j] += correlated_collision_integral( ( ztallymin - prevpos.z ) / fabs( dir.z ), step );
                    }
                }
            }
        }

        printDescriptor();
        std::cout << "Analogmu_t c: " << c << 
                     " maxz: " << maxz << 
                     " dz: " << dz << 
                     " du: " << du << 
                     " numsamples: " << numsamples << 
                     " numCollisionOrders: " << numCollisionOrders << 
                     " numMoments: " << numMoments << 
                     std::endl;

        std::cout << "CollisionDensity:\n";
        for( size_t i = 0; i < num_spatial_bins; ++i )
        {
            std::cout << double( collisionDensity[i] ) / ( double( numsamples ) * dz ) << " ";
        }
        std::cout << std::endl;

        std::cout << "Fluence:\n";
        for( size_t i = 0; i < num_spatial_bins; ++i )
        {
            std::cout << double( fluence[i] ) / ( double( numsamples ) * dz ) << " ";
        }
        std::cout << std::endl;

        delete [] collisionDensity;
        delete [] fluence;
    }
};
