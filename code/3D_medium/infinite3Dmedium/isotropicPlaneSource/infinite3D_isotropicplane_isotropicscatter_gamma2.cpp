#include "../GENinfinite3Dmedium.h"

class IsotropicGamma2InfiniteMedium : public GENInfiniteMedium
{
public:
    double m_mut;
    IsotropicGamma2InfiniteMedium(const double mut) : m_mut(mut){}

    double samplestep() { return Gamma2Step(); }
    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return isotropicDir();
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic plane source isotropic scattering Gamma2 random flight, mu_t = 1.0" << std::endl;
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

    sampler.isotropicPlaneSourceAnalogMut( c, maxz, dz, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}