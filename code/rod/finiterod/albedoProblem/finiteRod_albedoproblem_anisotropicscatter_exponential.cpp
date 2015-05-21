#include "../GENFiniteRod.h"

class HalfRodAlbedoProblemAnisotropicScatterExponential : public GENFiniteRod
{
public:
    double m_mu_t;
    double m_probF; // probability of scattering forward (F/c)
    double m_g;
    HalfRodAlbedoProblemAnisotropicScatterExponential(const double mu_t, const double in_probF, const double in_g ) : m_mu_t( mu_t ), m_probF( in_probF ), m_g( in_g ){}

    double init( double & pos, double & dir )
    {
        pos = 0.0;
        dir = 1.0;

        return 1.0;
    }
    double samplestep() { return -log( RandomReal() ) / m_mu_t; } // exponential sampling
    double sampledir( const double prevDir ) // Anisotropic scattering
    {
        if( RandomReal() < m_probF ) {return prevDir;} else {return -prevDir;};
    }

    void printDescriptor()
    {
        std::cout << "Finite rod, albedo problem, anisotropic scattering, exponential random flight, mu_t: " << m_mu_t << " g: " << m_g << std::endl;
    }
};

int main( int argc, char** argv )
{
    if( argc != 9 )
    {
        std::cout << "usage: finiteRodAlbedo c mu_t length dx numsamples numCollisionOrders numMoments g \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mu_t = StringToNumber<double>( std::string( argv[2] ) );
    double length = StringToNumber<double>( std::string( argv[3] ) );
    double dx = StringToNumber<double>( std::string( argv[4] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[5] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[6] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[7] ) );
    double g = StringToNumber<double>( std::string( argv[8] ) );

    const double F = 0.5 * ( 1.0 + g ) * c;

    HalfRodAlbedoProblemAnisotropicScatterExponential sampler( mu_t, F/c, g );

    sampler.RodAlbedoEstimatorAnalog( c, length, dx, numsamples, numCollisionOrders, numMoments );

    return 0;
}