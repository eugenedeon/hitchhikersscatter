#include "../GENinfinite3DmediumGRT.h"

// correlated emission - Cauchy scattering

class IsotropicCauchyCInfiniteMedium : public GENInfiniteMedium
{
public:
    
    IsotropicCauchyCInfiniteMedium() {}

    double sample_correlated_step() { return CauchyStep(); } 
    double sample_uncorrelated_step() { return CauchyStep(); } 
    double sigma_tc( const double s )
    {
        return 2.0 / ( ( 1.0 + s * s ) * ( M_PI - 2.0 * atan(s) ) );
    }
    double sigma_tu( const double s )
    {
        return 2.0 / ( ( 1.0 + s * s ) * ( M_PI - 2.0 * atan(s) ) );
    }
    double correlated_collision_integral( const double s1, const double s2 )
    {
        return Log((Pi - 2*atan(s1))/(Pi - 2*atan(s2)));
    }

    double uncorrelated_collision_integral( const double s1, const double s2 )
    {
        return Log((Pi - 2*atan(s1))/(Pi - 2*atan(s2)));
    }

    Vector3 sampledir( const Vector3& prevDir ) // rayleigh scattering
    {
        return RayleighDir(prevDir);
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic point source Rayleigh scattering Cauchy random flight.\n";
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

    IsotropicCauchyCInfiniteMedium sampler;

    sampler.isotropicPointSourceAnalogCollisionEstimator( c, maxr, dr, maxt, dt, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}