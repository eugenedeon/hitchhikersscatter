#include "../GENinfinite3DmediumGRT.h"

class IsoScatterGamma2CInfiniteMedium : public GENInfiniteMedium
{
public:
    double m_mfp;
    IsoScatterGamma2CInfiniteMedium( const double mfp ) : m_mfp(mfp){}

    double sample_correlated_step() { return Gamma2Step(); } 
    double sample_uncorrelated_step() { return Gamma2Step(); } 
    double sigma_tc( const double s )
    {
        return s / ( 1.0 + s );
    }
     double sigma_tu( const double s )
    {
        return s / ( 1.0 + s );
    }
    double correlated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + log(1.0 + s1) - log(1.0 + s2);
    }

    double uncorrelated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + log(1.0 + s1) - log(1.0 + s2);
    }

    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return isotropicDir();
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic point source Linearly-anisotropic scattering Gamma2C random flight, mfp = " << m_mfp << std::endl;
    }
};

int main( int argc, char** argv )
{
    srand48( time(NULL) );

    if( argc != 11 )
    {
        std::cout << "usage: infiniteMedium c mfp maxr dr du numsamples numCollisionOrders numMoments maxt dt \n";
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
    double maxt = StringToNumber<double>( std::string( argv[9] ) );
    double dt = StringToNumber<double>( std::string( argv[10] ) );

    IsoScatterGamma2CInfiniteMedium sampler( mfp );

    sampler.isotropicPointSourceAnalogCollisionEstimator( c, maxr, dr, maxt, dt, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}