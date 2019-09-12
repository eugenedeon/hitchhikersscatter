#include <Hfunction.h>
#include <fresnel.h>

double psi3D( const double u, const double c )
{
    return c * 0.5;
}

// [Williams 2006 - The albedo problem with Fresnel reflection ] - Eq.(64)
int main( int argc, char** argv )
{
    if( argc != 4 )
    {
        std::cout << "usage: deterministic c eta mui \n";
        exit( -1 );
    }

    const double c = StringToNumber<double>( std::string( argv[1] ) );
    const double eta = StringToNumber<double>( std::string( argv[2] ) );
    const double ui = StringToNumber<double>( std::string( argv[3] ) );

    const size_t QUADRATURE_ORDER( 100 );
    const double eps( 1e-8 );

    double I0[QUADRATURE_ORDER];
    double lastI0[QUADRATURE_ORDER];
    double approx0[QUADRATURE_ORDER];
    double Hs[QUADRATURE_ORDER];
    double Fs[QUADRATURE_ORDER];

    const double rui = fresnel::refractCosine( ui, 1.0, eta );
    const double Hrui = Hfunction100( rui, c, psi3D );

    // initial estimate:
    for( size_t i = 0; i < QUADRATURE_ORDER; ++i )
    {
        const double uo( Gauss100xs[i] );

        Hs[i] = Hfunction100( uo, c, psi3D );
        Fs[i] = fresnel::DielectricR( uo, 1.0 / eta );

        I0[i] = 0.5 * c * Hs[i] * Hrui * ui * ( 1.0 - fresnel::DielectricR( ui, eta ) ) / ( rui + uo );
        approx0[i] = I0[i];
    }

    for( size_t iter_i = 0; iter_i < 100; ++iter_i )
    {
        for( size_t i = 0; i < QUADRATURE_ORDER; ++i )
        {
            lastI0[i] = I0[i];
        }

        for( size_t i = 0; i < QUADRATURE_ORDER; ++i )
        {
            double sum( 0.0 );
            for( size_t j = 0; j < QUADRATURE_ORDER; ++j )
            {
                const double u( Gauss100xs[j] );
                sum += Gauss100ws[j] * u * Fs[j] * Hs[j] * lastI0[j] / ( Gauss100xs[i] + u );
            }
            I0[i] = approx0[i] + 0.5 * c * Hs[i] * sum;
        }
    }

    for( size_t i = 0; i < QUADRATURE_ORDER; ++i )
    {
        std::cout << I0[i] << " ";
    }
    std::cout << std::endl;
}