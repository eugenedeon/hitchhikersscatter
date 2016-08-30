#include "../GENHalfspace.h"

class IsotropicExponentialHalfspaceWhiteSkyAlbedo : public GENHalfspace
{
public:
    double mu_t;
    IsotropicExponentialHalfspaceWhiteSkyAlbedo(const double mu_t) : mu_t( mu_t ){}

    Vector3 initPos() { return Vector3( 0.0, 0.0, 0.0); }
    Vector3 initDir() { return lambertDir(); }

    double samplestep() { return -log( RandomReal() ) / mu_t; } // exponential sampling
    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return isotropicDir();
    }

    void printDescriptor()
    {
        std::cout << "Half-space white-sky albedo problem, isotropic scattering, exponential random flight, mu_t = " << mu_t << std::endl;
    }
};

int main( int argc, char** argv )
{
    if( argc != 9 )
    {
        std::cout << "usage: whiteskyalbedo c mu_t du dz maxz numsamples numCollisionOrders numMoments \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mu_t = StringToNumber<double>( std::string( argv[2] ) );
    double du = StringToNumber<double>( std::string( argv[3] ) );
    double dz = StringToNumber<double>( std::string( argv[4] ) );
    double maxz = StringToNumber<double>( std::string( argv[5] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[6] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[7] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[8] ) );

    IsotropicExponentialHalfspaceWhiteSkyAlbedo sampler( mu_t );

    sampler.HalfSpaceEstimatorAnalog( c, du, dz, maxz, numsamples, numCollisionOrders, numMoments );

    return 0;
}