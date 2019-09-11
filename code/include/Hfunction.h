#pragma once

#include <util.h>

// coordinate z > 0, absorption level c, characteristic function Psi( mu, c )
// evaluated iteratively using a 100-point Gauss quadrature
double Hfunction100( const double in_z, const double in_c, double (*in_psi)( const double, const double ) )
{
    const size_t QUADRATURE_ORDER( 100 );
    const double eps( 1e-8 );

    double Hm[QUADRATURE_ORDER];
    double lastHm[QUADRATURE_ORDER];
    double Psis[QUADRATURE_ORDER];

    for( size_t i = 0; i < QUADRATURE_ORDER; ++i )
    {
        Hm[i] = 1.0;
        Psis[i] = in_psi( Gauss100xs[i], in_c );
    }

    double result( 1.0 );
    double last_result( 0.0 );
    while( fabs( last_result - result ) > eps )
    {
        last_result = result;

        for( size_t i = 0; i < QUADRATURE_ORDER; ++i )
        {
            lastHm[i] = Hm[i];
        }

        for( size_t i = 0; i < QUADRATURE_ORDER; ++i )
        {
            double sum( 0.0 );
            for( size_t j = 0; j < QUADRATURE_ORDER; ++j )
            {
                sum += Gauss100ws[j] * lastHm[j] * Psis[j] / ( Gauss100xs[i] + Gauss100xs[j] );
            }  
            Hm[i] = 1.0 / ( 1.0 - Gauss100xs[i] * sum );
        }

        double sum( 0.0 );
        for( size_t j = 0; j < QUADRATURE_ORDER; ++j )
        {
            sum += Gauss100ws[j] * lastHm[j] * Psis[j] / ( in_z + Gauss100xs[j] );
        }  
        result = 1.0 / ( 1.0 - in_z * sum );
    }

    return result;
}