#pragma once

#include <vector.h>

//////////////////////////////////////////////////////////////////////////////////
// Random numbers
//////////////////////////////////////////////////////////////////////////////////

double RandomReal()
{
    return drand48();
}

double RandomReal( double a, double b )
{
    double r = drand48();
    return r * b + ( 1.0 - r ) * a;
}

// TODO: Box Muller can be numerically unstable - check bookmarked article
double RandomGauss()
{
  return sqrt(2.0)*cos(2*M_PI*RandomReal())*sqrt(-log(RandomReal()));
}

// normalize a 1D Gaussian distribution over [0,infty] such that the mean is 1.0
double GaussStep()
{
  return fabs( RandomGauss() ) / sqrt( 2.0 / M_PI );
}

double PearsonStep()
{
  return 1.0;
}

double BoxStep()
{
  return RandomReal(0.0, 2.0);
}

double Chik2Step()
{
  return (2.0 * sqrt( log( 1.0 / ( 1.0 - RandomReal() ) ) ) ) / sqrt(M_PI);
}

// Dagum distribution p = 2, a = 2
double Dagum22step()
{
  const float x = RandomReal();
  return (4*Sqrt(-(Sqrt(x)/(Power(M_PI,2)*(-1 + x))) - 
       x/(Power(M_PI,2)*(-1 + x))))/3.0;
}

double LogCauchyStep( const double a )
{
  return exp(-( a * 1.0 / tan( M_PI * RandomReal() ) ) );
}

//////////////////////////////////////////////////////////////////////////////////
// Random Directions
//////////////////////////////////////////////////////////////////////////////////

Vector3 diskSample( const double radius)
{
  const double phi = RandomReal( 0.0, 2 * M_PI );
  return radius * sqrt( RandomReal() ) * Vector3( cos( phi ), sin( phi ), 0.0 );
}

Vector3 isotropicDir()
{
    const double w = RandomReal( -1.0, 1.0 );
    const double p = RandomReal( 0.0, 2.0 * M_PI );
    const double s = sqrt( 1.0 - w * w );
    return Vector3( w, s * cos(p), s * sin(p) );
}

Vector2 isotropicDirFlatland()
{
    const double p = RandomReal( 0.0, 2.0 * M_PI );
    return Vector2( cos(p), sin(p) );
}

Vector4 isotropicDir4D()
{
  const double phi1 = 2 * M_PI * RandomReal();
  const double phi2 = 2 * M_PI * RandomReal();
  const double e1 = sqrt( -log( RandomReal() ) );
  const double e2 = sqrt( -log( RandomReal() ) );
  const double g1 = sin( phi1 ) * e1;
  const double g2 = cos( phi1 ) * e1;
  const double g3 = sin( phi2 ) * e2;
  const double g4 = cos( phi2 ) * e2;
  return Normalize( Vector4( g1, g2, g3, g4 ) );
}

Vector3 lambertDir()
{
    const double w = sqrt( RandomReal() );
    const double p = RandomReal( 0.0, 2.0 * M_PI );
    const double s = sqrt( 1.0 - w * w );
    return Vector3( s * cos(p), s * sin(p), w ); // z axis is the normal
}

// Henyey-Greenstein sampling
// in: is the direction the particle is moving before scatter
//  (not the direction the particle comes FROM)
Vector3 hgDir( const Vector3& in, const double g )
{
    if( fabs( g ) < 1e-5 )
    {
      return isotropicDir();
    }

    double t = ( 1.0 - g * g ) / ( 1.0 - g + 2.0 * g * RandomReal() );
    double w = ( 1.0 + g * g - t * t ) / ( 2.0 * g );
    const double p = RandomReal( 0.0, 2.0 * M_PI );
    const double s = sqrt( 1.0 - w * w );
    Vector3 aligned = Vector3( w, s * cos(p), s * sin(p) );

    if( in != Vector3( 0.0, 1.0, 0.0 ) )
    {
      const Vector3 v2 = Normalize( Vector3( 0.0, 1.0, 0.0 ) - in * in.y );
      const Vector3 v3 = Cross( in, v2 );

      return aligned.x * in + aligned.y * v2 + aligned.z * v3;
    }
    else
    {
      const Vector3 v2 = Vector3( 1.0, 0.0, 0.0 );
      const Vector3 v3 = Vector3( 0.0, 0.0, 1.0 );

      return aligned.x * in + aligned.y * v2 + aligned.z * v3;
    }
}

// linear-anisotropic scattering
Vector3 linanisoDir( const Vector3& in, const double b )
{
    double w = ( -1.0 + sqrt( 1.0 - 2.0*b + b*b + 4.0*b*RandomReal() ) ) / b;
    const double p = RandomReal( 0.0, 2.0 * M_PI );
    const double s = sqrt( 1.0 - w * w );
    Vector3 aligned = Vector3( w, s * cos(p), s * sin(p) );

    if( in != Vector3( 0.0, 1.0, 0.0 ) )
    {
      const Vector3 v2 = Normalize( Vector3( 0.0, 1.0, 0.0 ) - in * in.y );
      const Vector3 v3 = Cross( in, v2 );

      return aligned.x * in + aligned.y * v2 + aligned.z * v3;
    }
    else
    {
      const Vector3 v2 = Vector3( 1.0, 0.0, 0.0 );
      const Vector3 v3 = Vector3( 0.0, 0.0, 1.0 );

      return aligned.x * in + aligned.y * v2 + aligned.z * v3;
    }
}

// forward-scattering VMF (spherical gaussian) scattering kernel sampling
// fails for k > 500 or so
Vector3 forward_vmfDir( const Vector3& in, const double k )
{
    if( fabs( k ) < 1e-5 ) // todo update this
    {
      return isotropicDir();
    }

    const double e = RandomReal();
    const double w = log( e * exp( k ) + ( 1.0 - e ) * exp( -k ) ) / k;
    const double p = RandomReal( 0.0, 2.0 * M_PI );
    const double s = sqrt( 1.0 - w * w );
    Vector3 aligned = Vector3( w, s * cos(p), s * sin(p) );

    if( in != Vector3( 0.0, 1.0, 0.0 ) )
    {
      const Vector3 v2 = Normalize( Vector3( 0.0, 1.0, 0.0 ) - in * in.y );
      const Vector3 v3 = Cross( in, v2 );

      return aligned.x * in + aligned.y * v2 + aligned.z * v3;
    }
    else
    {
      const Vector3 v2 = Vector3( 1.0, 0.0, 0.0 );
      const Vector3 v3 = Vector3( 0.0, 0.0, 1.0 );

      return aligned.x * in + aligned.y * v2 + aligned.z * v3;
    }
}

// backward-scattering VMF (spherical gaussian) scattering kernel sampling
// fails for k > 500 or so
Vector3 backward_vmfDir( const Vector3& in, const double k )
{
    if( fabs( k ) < 1e-5 ) // todo update this
    {
      return isotropicDir();
    }

    const double e = RandomReal();
    const double w = log( e * exp( k ) + ( 1.0 - e ) * exp( -k ) ) / k;
    const double p = RandomReal( 0.0, 2.0 * M_PI );
    const double s = sqrt( 1.0 - w * w );
    Vector3 aligned = Vector3( -w, s * cos(p), s * sin(p) );

    if( in != Vector3( 0.0, 1.0, 0.0 ) )
    {
      const Vector3 v2 = Normalize( Vector3( 0.0, 1.0, 0.0 ) - in * in.y );
      const Vector3 v3 = Cross( in, v2 );

      return aligned.x * in + aligned.y * v2 + aligned.z * v3;
    }
    else
    {
      const Vector3 v2 = Vector3( 1.0, 0.0, 0.0 );
      const Vector3 v3 = Vector3( 0.0, 0.0, 1.0 );

      return aligned.x * in + aligned.y * v2 + aligned.z * v3;
    }
}

// build a forward/backward combination of vMF phase functions such that the square of the mean
// cosine is g2, and the 
Vector3 vmfvmfDir( const Vector3& in, const double k_forward, const double k_backward, const double mix )
{
  if( RandomReal() < mix )
  {
    return backward_vmfDir( in, k_backward );
  }
  else
  {
    return forward_vmfDir( in, k_forward );
  }
}