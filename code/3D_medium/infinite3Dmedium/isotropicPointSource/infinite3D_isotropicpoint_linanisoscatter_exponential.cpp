#include "../GENinfinite3DmediumGRT.h"

class LinAnisoExponentialInfiniteMedium : public GENInfiniteMedium
{
public:
    double m_mfp;
    double m_b;
    LinAnisoExponentialInfiniteMedium( const double mfp, const double b ) : m_mfp(mfp), m_b(b){}

    double sample_correlated_step() { return -log( RandomReal() ) * m_mfp; }
    double sample_uncorrelated_step() { return -log( RandomReal() ) * m_mfp; }
    double sigma_tc( const double s ) { return 1 / m_mfp; }
    double sigma_tu( const double s ) { return 1 / m_mfp; }
    double correlated_collision_integral( const double s1, const double s2 ) { return ( s2 - s1) / m_mfp; }
    double uncorrelated_collision_integral( const double s1, const double s2 ) { return ( s2 - s1) / m_mfp; }

    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return linanisoDir( prevDir, m_b );
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic point source Linearly-anisotropic scattering exponential random flight, mfp = " << m_mfp << " b = " << m_b << std::endl;
    }
};

int main( int argc, char** argv )
{
    srand48( time(NULL) );

    if( argc != 10 )
    {
        std::cout << "usage: infiniteMedium c mfp maxr dr du numsamples numCollisionOrders numMoments b \n";
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
    double b = StringToNumber<double>( std::string( argv[9] ) );

    LinAnisoExponentialInfiniteMedium sampler( mfp, b );

    sampler.isotropicPointSourceAnalogCollisionEstimator( c, maxr, dr, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}