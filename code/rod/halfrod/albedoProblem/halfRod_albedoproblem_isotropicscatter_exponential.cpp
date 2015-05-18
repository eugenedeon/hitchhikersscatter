#include "../GENHalfRod.h"

class HalfRodAlbedoProblemIsotropicScatterExponential : public GENHalfRod
{
public:
    double m_mu_t;
    HalfRodAlbedoProblemIsotropicScatterExponential(const double mu_t) : m_mu_t( mu_t ){}

    double init( double & pos, double & dir )
    {
        pos = 0.0;
        dir = 1.0;

        return 1.0;
    }
    double samplestep() { return -log( RandomReal() ) / m_mu_t; } // exponential sampling
    double sampledir( const double prevDir ) // isotropic scattering
    {
        if( RandomReal() > 0.5 ) {return 1.0;} else {return -1.0;};
    }

    void printDescriptor()
    {
        std::cout << "Half rod, albedo problem, isotropic scattering, exponential random flight, mu_t: " << m_mu_t << std::endl;
    }
};

int main( int argc, char** argv )
{
    if( argc != 8 )
    {
        std::cout << "usage: halfRodAlbedo c mu_t maxx dx numsamples numCollisionOrders numMoments \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mu_t = StringToNumber<double>( std::string( argv[2] ) );
    double maxx = StringToNumber<double>( std::string( argv[3] ) );
    double dx = StringToNumber<double>( std::string( argv[4] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[5] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[6] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[7] ) );

    HalfRodAlbedoProblemIsotropicScatterExponential sampler( mu_t );

    sampler.halfRodEstimatorAnalog( c, maxx, dx, numsamples, numCollisionOrders, numMoments );

    return 0;
}