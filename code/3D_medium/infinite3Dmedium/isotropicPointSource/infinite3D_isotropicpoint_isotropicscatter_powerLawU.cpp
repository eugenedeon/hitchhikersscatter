#include "../GENinfinite3DmediumGRT.h"

// correlated emission - power law extinction

class IsotropicPowerLawUInfiniteMedium : public GENInfiniteMedium
{
public:
    double m_mfp;
    double m_a;
    IsotropicPowerLawUInfiniteMedium(const double mfp, const double a) : m_mfp(mfp), m_a(a){}

    double sample_correlated_step() { return PowerLawCorrelatedStep( m_mfp, m_a ); } 
    double sample_uncorrelated_step() { return PowerLawUncorrelatedStep( m_mfp, m_a ); } 
    double correlated_transmittance( const double s )
    {
        return PowerLawCorrelatedTransmittance( s, m_a, m_mfp );
    }
    double uncorrelated_transmittance( const double s )
    {
        return PowerLawUncorrelatedTransmittance( s, m_a, m_mfp );
    }
    Vector3 sampledir( const Vector3& prevDir ) // isotropic scattering
    {
        return isotropicDir();
    }

    double sigma_tc(const double s)
    {
        return (1.0 + m_a)/(m_a*m_mfp + s);
    }

    double correlated_collision_integral( const double s1, const double s2 )
    {
        return (1.0 + m_a) * log((m_a*m_mfp + s2)/(m_a*m_mfp + s1));
    }

    double uncorrelated_collision_integral( const double s1, const double s2 )
    {
        return (m_a) * log((m_a*m_mfp + s2)/(m_a*m_mfp + s1));
    }

    double sigma_tu(const double s)
    {
        return (m_a)/(m_a*m_mfp + s);
    }

    void printDescriptor()
    {
        std::cout << "Infinite 3D isotropic point source isotropic scattering PowerLawU random flight, mfp = " << m_mfp << " a = " << m_a << std::endl;
    }
};

int main( int argc, char** argv )
{
    srand48( time(NULL) );

    if( argc != 10 )
    {
        std::cout << "usage: infiniteMedium c mfp maxr dr du numsamples numCollisionOrders numMoments a \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double mfp = StringToNumber<double>( std::string( argv[2] ) );
    double maxr = StringToNumber<double>( std::string( argv[3] ) );
    double dr = StringToNumber<double>( std::string( argv[4] ) );
    double du = StringToNumber<double>( std::string( argv[5] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[6] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[7] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[8] ) );
    double a = StringToNumber<double>( std::string( argv[9] ) );

    IsotropicPowerLawUInfiniteMedium sampler(mfp, a);

    sampler.isotropicPointSourceAnalogCollision( c, maxr, dr, du, numsamples, numCollisionOrders, numMoments );

    return 0;
}