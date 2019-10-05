#include "../GENinfinite3DmediumGRT.h"

class IsotropicGamma2InfiniteMedium : public GENInfiniteMedium
{
public:
    double m_mut;
    IsotropicGamma2InfiniteMedium(const double mut) : m_mut(mut){}

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

    double correlated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + log(1.0 + s1) - log(1.0 + s2);
    }
    double uncorrelated_collision_integral( const double s1, const double s2 )
    {
        return -s1 + s2 + log(1.0 + s1) - log(1.0 + s2);
    }

    double sigma_tc( const double s )
    {
        return s / ( 1.0 + s );
    }
    double sigma_tu( const double s )
    {
        return (1.0 + s) / ( 2.0 + s );
    }

    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return isotropicDir();
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic plane source isotropic scattering Gamma2 correlated random flight Track-length, mu_t = 1.0" << std::endl;
    }
};

int main( int argc, char** argv )
{
    srand48( time(NULL) );

    if( argc != 9 )
    {
        std::cout << "usage: infiniteMedium c mut maxz dz du numsamples numCollisionOrders numMoments \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mut = StringToNumber<double>( std::string( argv[2] ) );
    double maxz = StringToNumber<double>( std::string( argv[3] ) );
    double dz = StringToNumber<double>( std::string( argv[4] ) );
    double du = StringToNumber<double>( std::string( argv[5] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[6] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[7] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[8] ) );

    IsotropicGamma2InfiniteMedium sampler(mut);

    sampler.isotropicPlaneSourceAnalogTrackLength( c, maxz, dz, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}