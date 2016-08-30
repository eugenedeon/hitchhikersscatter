#include "../GENHalfspace.h"

class IsotropicExponentialHalfspaceDeltaAlbedo : public GENHalfspace
{
public:
    double mu_t;
    double m_ui;
    IsotropicExponentialHalfspaceDeltaAlbedo(const double mu_t, const double ui) : mu_t( mu_t ), m_ui( ui ){}

    Vector3 initPos() { return Vector3( 0.0, 0.0, 0.0); }
    Vector3 initDir() { return Vector3( 0.0, sqrt( 1.0 - m_ui * m_ui ), m_ui ); }

    double samplestep() { return -log( RandomReal() ) / mu_t; } // exponential sampling
    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return isotropicDir();
    }

    void printDescriptor()
    {
        std::cout << "Half-space delta albedo problem, isotropic scattering, exponential random flight.  ui: " << m_ui << " mu_t = " << mu_t << std::endl;
    }
};

int main( int argc, char** argv )
{
    if( argc != 10 )
    {
        std::cout << "usage: albedo_delta c mu_t ui du dz maxz numsamples numCollisionOrders numMoments \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mu_t = StringToNumber<double>( std::string( argv[2] ) );
    double ui = StringToNumber<double>( std::string( argv[3] ) );
    double du = StringToNumber<double>( std::string( argv[4] ) );
    double dz = StringToNumber<double>( std::string( argv[5] ) );
    double maxz = StringToNumber<double>( std::string( argv[6] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[7] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[8] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[9] ) );

    IsotropicExponentialHalfspaceDeltaAlbedo sampler( mu_t, ui );

    sampler.HalfSpaceEstimatorAnalog( c, du, dz, maxz, numsamples, numCollisionOrders, numMoments );

    return 0;
}