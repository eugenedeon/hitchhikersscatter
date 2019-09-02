#pragma once

#include <util.h>
#include <vector.h>

// exact [Dunkle 1963]
double SmoothDielectricHemisphericalAlbedo( const double n )
{
    if( n <= 0.0 )
    {
        return 0.0;
    }
    if( n < 1.0 )
    {
        return 1.0 - n * n * ( 1.0 - SmoothDielectricHemisphericalAlbedo( 1.0 / n ) );
    }
    if( n < 1.00000000001 )
    {
        return 0.0;
    }
    if( n < 1.001)
    {
        return (-1 + n)/3. - (79*Power(-1 + n,3))/60. + 
            (Power(-1 + n,4)*(37 + 100*Log(2) - 100*Log(-1 + n)))/160. + 
            (Power(-1 + n,2)*(19 - 12*Log(2) + 12*Log(-1 + n)))/24.;
    }
    return 0.5 + ((-1 + n)*(1 + 3*n))/(6.*Power(1 + n,2)) - 
        (2*Power(n,3)*(-1 + 2*n + Power(n,2)))/((1 + Power(n,2))*(-1 + Power(n,4))) + 
        (8*Power(n,4)*(1 + Power(n,4))*Log(n))/((1 + Power(n,2))*Power(-1 + Power(n,4),2)) + 
        (Power(n,2)*Power(-1 + Power(n,2),2)*Log((-1 + n)/(1 + n)))/Power(1 + Power(n,2),3);
}