#include "../GENinfinite3DmediumGRT.h"

class IsotropicExponentialInfiniteMedium : public GENInfiniteMedium
{
public:
    double m_mfp;
    IsotropicExponentialInfiniteMedium(const double mfp) : m_mfp(mfp){}

    double sample_correlated_step() { return -log( RandomReal() ) * m_mfp; }
    double sample_uncorrelated_step() { return -log( RandomReal() ) * m_mfp; }
    double sigma_tc( const double s ) { return 1 / m_mfp; }
    double sigma_tu( const double s ) { return 1 / m_mfp; }
    double correlated_collision_integral( const double s1, const double s2 ) { return ( s2 - s1) / m_mfp; }
    double uncorrelated_collision_integral( const double s1, const double s2 ) { return ( s2 - s1) / m_mfp; }
    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return lambertSphereDir(prevDir);
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic point source LambertSphere scattering exponential random flight, mfp = " << m_mfp << std::endl;
    }
};

int main( int argc, char** argv )
{
    srand48( time(NULL) );

    if( argc != 9 )
    {
        std::cout << "usage: infiniteMedium c mfp maxr dr du numsamples numCollisionOrders numMoments \n";
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

    IsotropicExponentialInfiniteMedium sampler(mfp);

    sampler.isotropicPointSourceAnalogCollisionEstimator( c, maxr, dr, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}