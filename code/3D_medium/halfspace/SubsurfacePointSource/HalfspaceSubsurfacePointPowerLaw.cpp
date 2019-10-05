#include "../GENHalfspaceGRT.h"

class HalfspaceSubsurfacePointPowerLaw : public GENHalfspace
{
public:
    double m_mfp;
    double m_depth;
    double m_a;
    HalfspaceSubsurfacePointPowerLaw( const double mfp, const double a, const double depth ) : m_mfp( mfp ), m_depth( depth ), m_a( a ) {}

    Vector3 initPos() { return Vector3( 0.0, 0.0, m_depth ); }
    Vector3 initDir() { return isotropicDir(); }

    double sample_correlated_step() { return PowerLawCorrelatedStep( m_mfp, m_a ); } 
    double sample_uncorrelated_step() { return PowerLawUncorrelatedStep( m_mfp, m_a ); } 
    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return isotropicDir();
    }

    void printDescriptor()
    {
        std::cout << "Half-space subsurface isotropic point, isotropic scattering, PowerLaw random flight.  depth: " << m_depth << " mfp = " << m_mfp << " a = " << m_a << std::endl;
    }
};

int main( int argc, char** argv )
{
    if( argc != 10 )
    {
        std::cout << "usage: HalfspaceSubsurfacePointPowerLaw c mfp depth dr maxr dz maxz numsamples a \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mfp = StringToNumber<double>( std::string( argv[2] ) );
    double depth = StringToNumber<double>( std::string( argv[3] ) );
    double dr = StringToNumber<double>( std::string( argv[4] ) );
    double maxr = StringToNumber<double>( std::string( argv[5] ) );
    double dz = StringToNumber<double>( std::string( argv[6] ) );
    double maxz = StringToNumber<double>( std::string( argv[7] ) );
    double a = StringToNumber<double>( std::string( argv[9] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[8] ) );

    HalfspaceSubsurfacePointPowerLaw sampler( mfp, a, depth );

    sampler.HalfSpace2DEstimatorAnalog( c, dz, maxz, dr, maxr, numsamples );

    return 0;
}