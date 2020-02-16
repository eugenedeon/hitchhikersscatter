#include "../GENinfinite3DmediumGRT.h"

// correlated emission - Gamma-6 scattering

class IsotropicGamma6CInfiniteMedium : public GENInfiniteMedium
{
public:
    
    IsotropicGamma6CInfiniteMedium() {}

    double sample_correlated_step() { return -log( RandomReal() * RandomReal() * RandomReal() * RandomReal() * RandomReal() * RandomReal() ); } 
    double sample_uncorrelated_step() { return -log( RandomReal() * RandomReal() * RandomReal() * RandomReal() * RandomReal() * RandomReal() ); } 
    double sigma_tc( const double s )
    {
        return Power(s,5)/(120 + s*(120 + s*(60 + s*(20 + s*(5 + s)))));
    }
     double sigma_tu( const double s )
    {
        return Power(s,5)/(120 + s*(120 + s*(60 + s*(20 + s*(5 + s)))));
    }
    double correlated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + Log(120 + s1*(120 + s1*(60 + s1*(20 + s1*(5 + s1))))) - 
   Log(120 + s2*(120 + s2*(60 + s2*(20 + s2*(5 + s2)))));
    }

    double uncorrelated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + Log(120 + s1*(120 + s1*(60 + s1*(20 + s1*(5 + s1))))) - 
   Log(120 + s2*(120 + s2*(60 + s2*(20 + s2*(5 + s2)))));
    }

    Vector3 sampledir( const Vector3& prevDir ) // rayleigh scattering
    {
        return RayleighDir(prevDir);
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic point source Rayleigh scattering Gamma-6 random flight.\n";
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

    IsotropicGamma6CInfiniteMedium sampler;

    sampler.isotropicPointSourceAnalogCollisionEstimator( c, maxr, dr, maxt, dt, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}