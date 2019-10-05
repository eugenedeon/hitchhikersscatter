#include "../GENhalfSpaceFlatland.h"

class IsotropicExponentialHalfSpace : public GENHalfSpace
{
public:
    double m_mut;
    IsotropicExponentialHalfSpace(const double mut) : m_mut(mut){}

    double samplestep() { return -log( RandomReal() ) / m_mut; }
    Vector2 sampledir( const Vector2& prevDir ) // isotropic scattering
    {
        return isotropicDirFlatland();
    }

    void printDescriptor()
    {
        std::cout << "Half-space flatland normally-incidence searchlight isotropic scattering exponential random flight, mu_t = " << m_mut << std::endl;
    }
};

int main( int argc, char** argv )
{
    srand48( time(NULL) );

    if( argc != 7 )
    {
        std::cout << "usage: halfSpace c mut maxr dr du numsamples \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mut = StringToNumber<double>( std::string( argv[2] ) );
    double maxr = StringToNumber<double>( std::string( argv[3] ) );
    double dr = StringToNumber<double>( std::string( argv[4] ) );
    double du = StringToNumber<double>( std::string( argv[5] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[6] ) );

    IsotropicExponentialHalfSpace sampler(mut);

    sampler.searchLightNormalIncidenceAnalog( c, maxr, dr, du, numsamples );

    return 0;
}