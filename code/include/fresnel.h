#pragma once

#include <util.h>
#include <vector.h>

//////////////////////////////////////////////////////////////////////////////////
// Fresnel
//////////////////////////////////////////////////////////////////////////////////

namespace fresnel
{
  double evalF( const double g, const double c )
  {
    return (pow(-c + g,2)*(1 + pow(-1 + c*(c + g),2)/pow(1 + c*(-c + g),2)))/
      (2.*pow(c + g,2));
  }

  // costheta = 1 is normal incidence
  // eta is ratio of new medium / current medium
  double DielectricR( const double costheta, const double eta )
  {
    const double sqrtinput = eta * eta - 1.0 + costheta * costheta;
    if( sqrtinput <= 0.0 )
    {
      return 1.0;
    }
    else
    {
      return evalF( sqrt( std::max( 0.0, sqrtinput ) ), costheta );
    }
  }

  Vector3 refract( const Vector3 in, const Vector3 normal, const double etai, const double etat )
  {
    return (etai*(in - normal*Dot(in,normal)))/etat - 
      normal*Sqrt(1 - (pow(etai,2)*(1 - pow(Dot(in,normal),2)))/pow(etat,2));
  }

  double refractCosine( const double ui, const double etai, const double etao )
  {
    return sqrt( 1.0 - etai * etai * ( 1.0 - ui * ui ) / ( etao * etao ) );
  }

  double ConductorR( const double costheta, const double etai, const double eta, const double k )
  {
    return ((Power(etai,2) - 2*Power(costheta,2)*Power(etai,2) + Power(costheta,4)*Power(etai,2) + 
       Power(costheta,2)*Sqrt(4*Power(eta,2)*Power(k,2) + 
          Power(Power(eta,2) + (-1 + Power(costheta,2))*Power(etai,2) - Power(k,2),2)))*
     (Power(costheta,2)*Power(etai,2) + Sqrt(4*Power(eta,2)*Power(k,2) + 
         Power(Power(eta,2) + (-1 + Power(costheta,2))*Power(etai,2) - Power(k,2),2)) - 
       Sqrt(2)*costheta*etai*Sqrt(Power(eta,2) - Power(etai,2) + Power(costheta,2)*Power(etai,2) - 
          Power(k,2) + Sqrt(4*Power(eta,2)*Power(k,2) + 
            Power(Power(eta,2) + (-1 + Power(costheta,2))*Power(etai,2) - Power(k,2),2)))))/
   ((Power(costheta,2)*Power(etai,2) + Sqrt(4*Power(eta,2)*Power(k,2) + 
         Power(Power(eta,2) + (-1 + Power(costheta,2))*Power(etai,2) - Power(k,2),2)) + 
       Sqrt(2)*costheta*etai*Sqrt(Power(eta,2) - Power(etai,2) + Power(costheta,2)*Power(etai,2) - 
          Power(k,2) + Sqrt(4*Power(eta,2)*Power(k,2) + 
            Power(Power(eta,2) + (-1 + Power(costheta,2))*Power(etai,2) - Power(k,2),2))))*
     (Power(-1 + Power(costheta,2),2)*Power(etai,2) + 
       Power(costheta,2)*Sqrt(Power(eta,4) + 
          Power(-((-1 + Power(costheta,2))*Power(etai,2)) + Power(k,2),2) + 
          2*Power(eta,2)*((-1 + Power(costheta,2))*Power(etai,2) + Power(k,2))) - 
       Sqrt(2)*costheta*(-1 + Power(costheta,2))*etai*
        Sqrt(Power(eta,2) - Power(etai,2) + Power(costheta,2)*Power(etai,2) - Power(k,2) + 
          Sqrt(Power(eta,4) + Power(-((-1 + Power(costheta,2))*Power(etai,2)) + Power(k,2),2) + 
            2*Power(eta,2)*((-1 + Power(costheta,2))*Power(etai,2) + Power(k,2))))));
  }
}

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