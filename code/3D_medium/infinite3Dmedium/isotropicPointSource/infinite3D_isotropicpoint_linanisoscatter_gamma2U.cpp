#include "../GENinfinite3DmediumGRT.h"

class LinAnisoGamma2UInfiniteMedium : public GENInfiniteMedium
{
public:
    double m_mfp;
    double m_b;
    LinAnisoGamma2UInfiniteMedium( const double mfp, const double b ) : m_mfp(mfp), m_b(b){}

    double sample_correlated_step() { return Gamma2Step(); }
    double sample_uncorrelated_step() { return Gamma2ExtinctionStep(); }
    double correlated_transmittance( const double s )
    {
        return exp(-s) * ( 1.0 + s );
    }
    double uncorrelated_transmittance( const double s )
    {
        return 0.5 * exp(-s) * ( 2.0 + s );
    }

    double sigma_tc( const double s )
    {
        return s / ( 1.0 + s );
    }
    double sigma_tu( const double s )
    {
        return (1.0 + s) / ( 2.0 + s );
    }

    double correlated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + log(1.0 + s1) - log(1.0 + s2);
    }

    double uncorrelated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + log(2.0 + s1) - log(2.0 + s2);
    }

    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return linanisoDir( prevDir, m_b );
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic point source Linearly-anisotropic scattering Gamma2U random flight, mfp = " << m_mfp << " b = " << m_b << std::endl;
    }
};

int main( int argc, char** argv )
{
    srand48( time(NULL) );

    if( argc != 12 )
    {
        std::cout << "usage: infiniteMedium c mfp maxr dr du numsamples numCollisionOrders numMoments b maxt dt \n";
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
    double b = StringToNumber<double>( std::string( argv[11] ) );

    LinAnisoGamma2UInfiniteMedium sampler( mfp, b );

    sampler.isotropicPointSourceAnalogCollisionEstimator( c, maxr, dr, maxt, dt, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}