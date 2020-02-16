#include "../GENinfinite3DmediumGRT.h"

class LinAnisoGamma3CInfiniteMedium : public GENInfiniteMedium
{
public:
    double m_mfp;
    double m_b;
    LinAnisoGamma3CInfiniteMedium( const double mfp, const double b ) : m_mfp(mfp), m_b(b){}

    double sample_correlated_step() { return -log( RandomReal() * RandomReal() * RandomReal() ); } 
    double sample_uncorrelated_step() { return -log( RandomReal() * RandomReal() * RandomReal() ); } 
    double sigma_tc( const double s )
    {
        return s * s/(2 + s*(2 + s));
    }
     double sigma_tu( const double s )
    {
        return s * s/(2 + s*(2 + s));
    }
    double correlated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + Log(2 + 2*s1 + Power(s1,2)) - Log(2 + 2*s2 + Power(s2,2));
    }

    double uncorrelated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + Log(2 + 2*s1 + Power(s1,2)) - Log(2 + 2*s2 + Power(s2,2));
    }

    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return linanisoDir( prevDir, m_b );
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic point source Linearly-anisotropic scattering Gamma3C random flight, mfp = " << m_mfp << " b = " << m_b << std::endl;
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

    LinAnisoGamma3CInfiniteMedium sampler( mfp, b );

    sampler.isotropicPointSourceAnalogCollisionEstimator( c, maxr, dr, maxt, dt, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}