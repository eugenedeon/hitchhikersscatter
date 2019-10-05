#include "../GENinfinite3DmediumGRT.h"

// correlated emission - Gamma-2 scattering - HG phase

class IsotropicGamma2CInfiniteMedium : public GENInfiniteMedium
{
public:
    
    double m_g;
    IsotropicGamma2CInfiniteMedium(const double in_g) : m_g(in_g) {}

    double sample_correlated_step() { return Gamma2Step(); } 
    double sample_uncorrelated_step() { return Gamma2Step(); } 
    double correlated_transmittance( const double s )
    {
        return exp(-s) * ( 1.0 + s );
    }
    double uncorrelated_transmittance( const double s )
    {
        return exp(-s) * ( 1.0 + s );
    }
    double sigma_tc( const double s )
    {
        return s / ( 1.0 + s );
    }
     double sigma_tu( const double s )
    {
        return s / ( 1.0 + s );
    }
    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return hgDir(prevDir,m_g);
    }

    double correlated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + log(1.0 + s1) - log(1.0 + s2);
    }

    double uncorrelated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + log(1.0 + s1) - log(1.0 + s2);
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic point source HG g=" << m_g <<" Gamma-2 random flight.\n";
    }
};

int main( int argc, char** argv )
{
    srand48( time(NULL) );

    if( argc != 10 )
    {
        std::cout << "usage: infiniteMedium c mfp maxr dr du numsamples numCollisionOrders numMoments g \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mfp = StringToNumber<double>( std::string( argv[2] ) );
    double maxr = StringToNumber<double>( std::string( argv[3] ) );
    double dr = StringToNumber<double>( std::string( argv[4] ) );
    double du = StringToNumber<double>( std::string( argv[5] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[6] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[7] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[8] ) );
    double g = StringToNumber<double>( std::string( argv[9] ) );

    IsotropicGamma2CInfiniteMedium sampler(g);

    sampler.isotropicPointSourceAnalogCollision( c, maxr, dr, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}