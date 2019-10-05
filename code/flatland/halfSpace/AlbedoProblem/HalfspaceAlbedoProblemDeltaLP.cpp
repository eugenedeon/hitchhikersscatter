#include "../GENHalfspaceFlatlandGRT.h"

class IsotropicLPHalfspaceDeltaAlbedo : public GENHalfspace
{
public:
    double m_ui;

    double m_rm;
    double m_rp;
    double m_wp;

    double m_pa;
    double m_sigratio;
    double m_lambsum;

    // pa - volume fraction of phase A
    // sigratio - SigmaA / SigmaB
    // lambsum - lambdaA + lambdaB
    IsotropicLPHalfspaceDeltaAlbedo(const double pa, const double sigratio, const double lambsum, const double ui): m_pa( pa ), m_sigratio( sigratio), m_lambsum( lambsum ), m_ui( ui )
    {
        const double pb = 1.0 - pa;
        const double sigA = sigratio/(1.0 - pa + sigratio*pa);
        const double sigB = sigA / sigratio;
        const double lambA = pa * lambsum;
        const double lambB = lambsum - lambA;
        const double beta = (sigA-sigB)*(sigA-sigB) * pa * pb;
        const double sigtilda = pb * sigA + pa * sigB + 1.0 / lambA + 1.0 / lambB;
        m_rp = (1.0 + Sqrt(4*beta + Power(1.0 - sigtilda,2)) + sigtilda)/2.;
        m_rm = (1.0 - Sqrt(4*beta + Power(1.0 - sigtilda,2)) + sigtilda)/2.;
        m_wp = 1.0 - (-m_rm + sigtilda)/(-m_rm + m_rp);
    }

    Vector2 initPos() { return Vector2( 0.0, 0.0); }
    Vector2 initDir() { return Vector2( sqrt( 1.0 - m_ui * m_ui ), m_ui ); }

    double sample_correlated_step()
    { 
        if( RandomReal() < m_rm * (1.0 - m_wp) )
        {
            return -log( RandomReal() ) / m_rm;
        }
        else
        {
            return -log( RandomReal() ) / m_rp;
        }
    } 
    double sample_uncorrelated_step()
    {
        if( RandomReal() < m_wp )
        {
            return -log( RandomReal() ) / m_rp;
        }
        else
        {
            return -log( RandomReal() ) / m_rm;
        }
    }

    double cmfp(){ return 1; }
    double umfp(){ return ( m_rp + m_rm * m_wp - m_rp * m_wp )/( m_rm * m_rp ); }

    double sigma_t_c( const double s )
    {
        return m_rp + (exp(m_rp*s)*m_rm*(m_rm - m_rp)*(-1 + m_wp))/(exp(m_rp*s)*m_rm*(-1 + m_wp) - exp(m_rm*s)*m_rp*m_wp);
    }

    double sigma_t_u( const double s )
    {
        return ((m_rm*(1 - m_wp))/exp(m_rm*s) + (m_rp*m_wp)/exp(m_rp*s))/((1 - m_wp)/exp(m_rm*s) + m_wp/exp(m_rp*s));
    }

    Vector2 sampledir( const Vector2& prevDir ) // isotropic scattering
    {
        return isotropicDirFlatland();
    }

    void printDescriptor()
    {
        std::cout << "Half-space delta albedo problem, isotropic scattering, LP random flight.  ui: " << m_ui 
        << " rm = " << m_rm 
        << " rp = " << m_rp 
        << " wp = " << m_wp 
        << " pa = " << m_pa 
        << " sigratio = " << m_sigratio
        << " lambsum = " << m_lambsum
        << std::endl;
    }
};

int main( int argc, char** argv )
{
    if( argc != 12 )
    {
        std::cout << "usage: albedo_delta c vola ui du dz maxz numsamples numCollisionOrders numMoments sigratio lambsum \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double pa = StringToNumber<double>( std::string( argv[2] ) );
    double ui = StringToNumber<double>( std::string( argv[3] ) );
    double du = StringToNumber<double>( std::string( argv[4] ) );
    double dz = StringToNumber<double>( std::string( argv[5] ) );
    double maxz = StringToNumber<double>( std::string( argv[6] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[7] ) );
    size_t numCollisionOrders = StringToNumber<size_t>( std::string( argv[8] ) );
    size_t numMoments = StringToNumber<size_t>( std::string( argv[9] ) );
    double sigratio = StringToNumber<double>( std::string( argv[10] ) );
    double lambsum = StringToNumber<double>( std::string( argv[11] ) );

    IsotropicLPHalfspaceDeltaAlbedo sampler( pa, sigratio, lambsum, ui );

    sampler.HalfSpaceEstimatorAnalog( c, du, dz, maxz, numsamples, numCollisionOrders, numMoments );

    return 0;
}