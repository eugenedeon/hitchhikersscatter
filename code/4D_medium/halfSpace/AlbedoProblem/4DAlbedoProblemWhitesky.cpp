#include "../GENHalfSpace4D.h"

class IsotropicExponential4DHalfspaceWhiteSkyAlbedo : public GENHalfspace
{
public:
    double mu_t;
    IsotropicExponential4DHalfspaceWhiteSkyAlbedo(const double mu_t) : mu_t( mu_t ){}

    Vector4 initPos() { return Vector4( 0.0, 0.0, 0.0, 0.0); }
    Vector4 initDir() { return lambertDir4D(); }

    double samplestep() { return -log( RandomReal() ) / mu_t; } // exponential sampling
    Vector4 sampledir( const Vector4& prevDir ) // isotropic scattering
    {
        return isotropicDir4D();
    }

    void printDescriptor()
    {
        std::cout << "4D Half-space white-sky albedo problem, isotropic scattering, exponential random flight, mu_t = " << mu_t << std::endl;
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

    IsotropicExponential4DHalfspaceWhiteSkyAlbedo sampler( mu_t );

    sampler.HalfSpaceEstimatorAnalog( c, du, dz, maxz, numsamples, numCollisionOrders, numMoments );

    return 0;
}