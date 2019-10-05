#include "../GENHalfRod.h"

class HalfRodPointSourceIsotropicScatterExponential : public GENHalfRod
{
public:
    double m_x0;
    HalfRodPointSourceIsotropicScatterExponential(const double in_x0) : m_x0( in_x0 ){}

    double init( double & pos, double & dir )
    {
        pos = m_x0;
        dir = RandomReal() > 0.5 ? 1.0 : -1.0;

        return 1.0;
    }
    double samplestep() { return -log( RandomReal() ); } // exponential sampling
    double sampledir( const double prevDir ) // isotropic scattering
    {
        if( RandomReal() > 0.5 ) {return 1.0;} else {return -1.0;};
    }

    void printDescriptor()
    {
        std::cout << "Half rod, point source, isotropic scattering, exponential random flight, x0: " << m_x0 << std::endl;
    }
};

int main( int argc, char** argv )
{
    if( argc != 8 )
    {
        std::cout << "usage: halfRodPointSource c x0 maxx dx numsamples numCollisionOrders numMoments \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double x0 = StringToNumber<double>( std::string( argv[2] ) );
    double maxx = StringToNumber<double>( std::string( argv[3] ) );
    double dx = StringToNumber<double>( std::string( argv[4] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[5] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[6] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[7] ) );

    HalfRodPointSourceIsotropicScatterExponential sampler( x0 );

    sampler.halfRodEstimatorAnalog( c, maxx, dx, numsamples, numCollisionOrders, numMoments );

    return 0;
}