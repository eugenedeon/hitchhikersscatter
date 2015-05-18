#include "../GENinfiniteRod.h"

class InfiniteRodIsotropicPointIsotropicScatterExponential : public GENInfiniteRod
{
public:
    double m_mu_t;
    InfiniteRodIsotropicPointIsotropicScatterExponential(const double mu_t) : m_mu_t( mu_t ){}

    double init( double & pos, double & dir )
    {
        pos = 0.0;
        if( RandomReal() < 0.5 )
        {
            dir = 1.0;
        }
        else
        {
            dir = -1.0;
        }

        return 1.0;
    }
    double samplestep() { return -log( RandomReal() ) / m_mu_t; } // exponential sampling
    double sampledir( const double prevDir ) // isotropic scattering
    {
        if( RandomReal() > 0.5 ) {return 1.0;} else {return -1.0;};
    }

    void printDescriptor()
    {
        std::cout << "Infinite rod, isotropic point source, isotropic scattering, exponential random flight, mu_t: " << m_mu_t << std::endl;
    }
};

int main( int argc, char** argv )
{
    if( argc != 8 )
    {
        std::cout << "usage: infiniteRod_isotropicpoint_isotropicscatter_exponential c mu_t maxx dx numsamples numCollisionOrders numMoments \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mu_t = StringToNumber<double>( std::string( argv[2] ) );
    double maxx = StringToNumber<double>( std::string( argv[3] ) );
    double dx = StringToNumber<double>( std::string( argv[4] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[5] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[6] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[7] ) );

    InfiniteRodIsotropicPointIsotropicScatterExponential sampler( mu_t );

    sampler.infiniteRodEstimatorAnalog( c, maxx, dx, numsamples, numCollisionOrders, numMoments );

    return 0;
}