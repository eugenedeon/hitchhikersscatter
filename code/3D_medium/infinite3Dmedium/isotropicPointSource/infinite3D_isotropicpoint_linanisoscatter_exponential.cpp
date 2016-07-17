#include "../GENinfinite3Dmedium.h"

class LinAnisoExponentialInfiniteMedium : public GENInfiniteMedium
{
public:
    double m_mut;
    double m_b;
    LinAnisoExponentialInfiniteMedium(const double mut, const double b) : m_mut(mut), m_b(b){}

    double samplestep() { return -log( RandomReal() ) / m_mut; }
    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return linanisoDir(prevDir,m_b);
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3d isotropic point source isotropic scattering exponential random flight, mu_t = " << m_mut << " b = " << m_b << std::endl;
    }
};

int main( int argc, char** argv )
{
    srand48( time(NULL) );

    if( argc != 10 )
    {
        std::cout << "usage: infiniteMedium c mut maxr dr du numsamples numCollisionOrders numMoments b \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mut = StringToNumber<double>( std::string( argv[2] ) );
    double maxr = StringToNumber<double>( std::string( argv[3] ) );
    double dr = StringToNumber<double>( std::string( argv[4] ) );
    double du = StringToNumber<double>( std::string( argv[5] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[6] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[7] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[8] ) );
    double b = StringToNumber<double>( std::string( argv[9] ) );

    LinAnisoExponentialInfiniteMedium sampler(mut, b);

    sampler.isotropicPointSourceAnalogMut( c, maxr, dr, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}