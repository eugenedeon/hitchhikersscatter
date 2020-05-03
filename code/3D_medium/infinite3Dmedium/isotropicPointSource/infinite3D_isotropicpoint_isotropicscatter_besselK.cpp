#include "../GENinfinite3DmediumGRT.h"

// correlated emission - BesselK scattering

class IsotropicBesselKCInfiniteMedium : public GENInfiniteMedium
{
public:
    
    IsotropicBesselKCInfiniteMedium() {}

    double sample_correlated_step() { return BesselK0Step(); } 
    double sample_uncorrelated_step() { return BesselK0Step(); } 
    double sigma_tc( const double s )
    {
        return 0.0;
    }
    double sigma_tu( const double s )
    {
        return 0.0;
    }
    double correlated_collision_integral( const double s1, const double s2 )
    {
        return 0.0;
    }

    double uncorrelated_collision_integral( const double s1, const double s2 )
    {
        return 0.0;
    }

    Vector3 sampledir( const Vector3& prevDir ) // rayleigh scattering
    {
        return isotropicDir();
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic point source Rayleigh scattering BesselK random flight.\n";
    }
};

int main( int argc, char** argv )
{
    srand48( time(NULL) );

    if( argc != 10 )
    {
        std::cout << "usage: infiniteMedium c maxr dr maxt dt du numsamples numCollisionOrders numMoments \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double maxr = StringToNumber<double>( std::string( argv[2] ) );
    double dr = StringToNumber<double>( std::string( argv[3] ) );
    double maxt = StringToNumber<double>( std::string( argv[4] ) );
    double dt = StringToNumber<double>( std::string( argv[5] ) );
    double du = StringToNumber<double>( std::string( argv[6] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[7] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[8] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[9] ) );

    IsotropicBesselKCInfiniteMedium sampler;

    sampler.isotropicPointSourceAnalogCollisionEstimator( c, maxr, dr, maxt, dt, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}