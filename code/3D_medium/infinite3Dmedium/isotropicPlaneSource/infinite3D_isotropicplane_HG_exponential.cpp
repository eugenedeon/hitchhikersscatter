#include "../GENinfinite3DmediumGRT.h"

// correlated emission - exponential media - HG phase

class IsotropicExponentialInfiniteMedium : public GENInfiniteMedium
{
public:
    
    double m_mfp;
    double m_g;
    IsotropicExponentialInfiniteMedium(const double in_mfp, const double in_g) : m_mfp( in_mfp), m_g(in_g) {}

    double sample_correlated_step() { return -m_mfp * log( RandomReal() ); } 
    double sample_uncorrelated_step() { return -m_mfp * log( RandomReal() ); } 
    double correlated_transmittance( const double s )
    {
        return exp(-s / m_mfp);
    }
    double uncorrelated_transmittance( const double s )
    {
        return exp(-s / m_mfp);
    }
    double sigma_tc( const double s )
    {
        return 1.0 / m_mfp;
    }
     double sigma_tu( const double s )
    {
        return 1.0 / m_mfp;
    }
    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return hgDir(prevDir,m_g);
    }

    double correlated_collision_integral( const double s1, const double s2 )
    {
        return (-s1 + s2) / m_mfp ;
    }

    double uncorrelated_collision_integral( const double s1, const double s2 )
    {
        return (-s1 + s2) / m_mfp ;
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic point source HG g=" << m_g <<" exponential random flight.\n";
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

    IsotropicExponentialInfiniteMedium sampler(mfp,g);

    sampler.isotropicPlaneSourceAnalogCollision( c, maxr, dr, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}