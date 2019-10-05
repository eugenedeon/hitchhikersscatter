#include "../GENHalfspace.h"

class HalfspaceSubsurfacePointExponential : public GENHalfspace
{
public:
    double m_mu_t;
    double m_depth;
    HalfspaceSubsurfacePointExponential( const double mu_t, const double depth ) : m_mu_t( mu_t ), m_depth( depth ){}

    Vector3 initPos() { return Vector3( 0.0, 0.0, m_depth ); }
    Vector3 initDir() { return isotropicDir(); }

    double samplestep() { return -log( RandomReal() ) / m_mu_t; } 
    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return isotropicDir();
    }

    void printDescriptor()
    {
        std::cout << "Half-space subsurface isotropic point, isotropic scattering, exponential random flight.  depth: " << m_depth << " mu_t = " << m_mu_t << std::endl;
    }
};

int main( int argc, char** argv )
{
    if( argc != 10 )
    {
        std::cout << "usage: HalfspaceSubsurfacePointExponential c mu_t depth du dz maxz numsamples numCollisionOrders numMoments \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mu_t = StringToNumber<double>( std::string( argv[2] ) );
    double depth = StringToNumber<double>( std::string( argv[3] ) );
    double du = StringToNumber<double>( std::string( argv[4] ) );
    double dz = StringToNumber<double>( std::string( argv[5] ) );
    double maxz = StringToNumber<double>( std::string( argv[6] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[7] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[8] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[9] ) );

    HalfspaceSubsurfacePointExponential sampler( mu_t, depth );

    sampler.HalfSpaceEstimatorAnalog( c, du, dz, maxz, numsamples, numCollisionOrders, numMoments );

    return 0;
}