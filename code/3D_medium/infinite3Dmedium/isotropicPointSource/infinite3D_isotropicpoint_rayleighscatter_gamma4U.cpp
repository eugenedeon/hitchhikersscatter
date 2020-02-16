#include "../GENinfinite3DmediumGRT.h"

// Uncorrelated emission - Gamma-4 scattering

class IsotropicGamma4UInfiniteMedium : public GENInfiniteMedium
{
public:
    
    IsotropicGamma4UInfiniteMedium() {}

    double sample_correlated_step() { return -log( RandomReal() * RandomReal() * RandomReal() * RandomReal() ); }
    double sample_uncorrelated_step() {
        if( RandomReal() < 0.5 )
        {
            if( RandomReal() < 0.5 )
            {
                return -log( RandomReal() );
            }
            else
            {
                return -log( RandomReal() * RandomReal() );
            }
        }
        else
        {
            if( RandomReal() < 0.5 )
            {
                return -log( RandomReal() * RandomReal() * RandomReal() );
            }
            else
            {
                return -log( RandomReal() * RandomReal() * RandomReal() * RandomReal() );
            }
        }
    }
    double sigma_tc( const double s )
    {
        return Power(s,3)/(6 + s*(6 + s*(3 + s)));
    }
     double sigma_tu( const double s )
    {
        return (6 + s*(6 + s*(3 + s)))/(24 + s*(18 + s*(6 + s)));
    }
    double correlated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + Log(6 + s1*(6 + s1*(3 + s1))) - Log(6 + s2*(6 + s2*(3 + s2)));
    }

    double uncorrelated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + Log(24 + s1*(18 + s1*(6 + s1))) - Log(24 + s2*(18 + s2*(6 + s2)));
    }

    Vector3 sampledir( const Vector3& prevDir ) // rayleigh scattering
    {
        return RayleighDir(prevDir);
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic point source Rayleigh scattering Gamma-4U random flight.\n";
    }
};

int main( int argc, char** argv )
{
    srand48( time(NULL) );

    if( argc != 11 )
    {
        std::cout << "usage: infiniteMedium c mfp maxr dr maxt dt du numsamples numCollisionOrders numMoments \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mfp = StringToNumber<double>( std::string( argv[2] ) );
    double maxr = StringToNumber<double>( std::string( argv[3] ) );
    double dr = StringToNumber<double>( std::string( argv[4] ) );
    double maxt = StringToNumber<double>( std::string( argv[5] ) );
    double dt = StringToNumber<double>( std::string( argv[6] ) );
    double du = StringToNumber<double>( std::string( argv[7] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[8] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[9] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[10] ) );

    IsotropicGamma4UInfiniteMedium sampler;

    sampler.isotropicPointSourceAnalogCollisionEstimator( c, maxr, dr, maxt, dt, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}