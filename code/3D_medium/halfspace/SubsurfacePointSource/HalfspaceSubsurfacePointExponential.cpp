#include "../GENHalfspaceGRT.h"

class HalfspaceSubsurfacePointExponential : public GENHalfspace
{
public:
    double m_mu_t;
    double m_depth;
    HalfspaceSubsurfacePointExponential( const double mu_t, const double depth ) : m_mu_t( mu_t ), m_depth( depth ){}

    Vector3 initPos() { return Vector3( 0.0, 0.0, m_depth ); }
    Vector3 initDir() { return isotropicDir(); }

    double sample_correlated_step() { return -log( RandomReal() ) / m_mu_t; } 
    double sample_uncorrelated_step() { return -log( RandomReal() ) / m_mu_t; } 
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
    if( argc != 9 )
    {
        std::cout << "usage: HalfspaceSubsurfacePointExponential c mu_t depth dr maxr dz maxz numsamples \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mu_t = StringToNumber<double>( std::string( argv[2] ) );
    double depth = StringToNumber<double>( std::string( argv[3] ) );
    double dr = StringToNumber<double>( std::string( argv[4] ) );
    double maxr = StringToNumber<double>( std::string( argv[5] ) );
    double dz = StringToNumber<double>( std::string( argv[6] ) );
    double maxz = StringToNumber<double>( std::string( argv[7] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[8] ) );

    HalfspaceSubsurfacePointExponential sampler( mu_t, depth );

    sampler.HalfSpace2DEstimatorAnalog( c, dz, maxz, dr, maxr, numsamples );

    return 0;
}