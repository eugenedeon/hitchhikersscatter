#pragma once

#include <vector.h>

//////////////////////////////////////////////////////////////////////////////////
// Random numbers
//////////////////////////////////////////////////////////////////////////////////

// uniform random real in [0,1]
double RandomReal()
{
    return drand48();
}

// uniform random real in [a,b]
double RandomReal( double a, double b )
{
    double r = RandomReal();
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

// step-distrubion p(s) = s e^(-s)
double Gamma2Step()
{
  return -log( RandomReal() * RandomReal() );
}

// step-distrubion p(s) = 0.5 * s e^(-s) + 0.5 * e^(-s)
double Gamma2ExtinctionStep()
{
  if( RandomReal() < 0.5 )
  {
    return Gamma2Step();
  }
  else
  {
    return -log( RandomReal() );
  }
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

double Chid3CorrelatedStep()
{
  return Norm( (2.0 / (M_PI ))*Vector3(RandomGauss(),RandomGauss(),RandomGauss() ) );
}

// a = 3 distribution - using modified definition in current book
double Chi3CorrelatedStep()
{
  // this isn't right - not sure why the 0.98 is required
  return 0.98 * 1.436696977001332493512655869021542042257 * Norm( Vector3( RandomGauss(), RandomGauss(), RandomGauss() ) );
}

double Chid3UncorrelatedStep()
{
  if( RandomReal() < 0.5 )
  {
    const double xi = RandomReal();
    return ( sqrt( 2.0 * M_PI * xi - M_PI * xi * xi) * fabs(RandomGauss()) )/( 2.0 * sqrt(2) );
  }
  else
  {
    const double xi = RandomReal();
    return (sqrt(M_PI)*sqrt(log(1.0/(1.0 - xi))))/2.0;
  }
}

// Dagum distribution p = 2, a = 2
double Dagum22step()
{
  const float x = RandomReal();
  return (4*Sqrt(-(Sqrt(x)/(Power(M_PI,2)*(-1 + x))) - 
       x/(Power(M_PI,2)*(-1 + x))))/3.0;
}

double CauchyStep()
{
  return tan( ( M_PI * RandomReal() ) / 2.0 );
}

double BesselK0Step()
{
  return -( cos( ( M_PI * RandomReal() ) / 2.0 ) * log( RandomReal() ) );
}

double LogCauchyStep( const double a )
{
  return exp(-( a * 1.0 / tan( M_PI * RandomReal() ) ) );
}

double PowerLawCorrelatedStep( const double mfp, const double a )
{
  return a * mfp * ( pow( RandomReal(), -1.0 / ( 1.0 + a ) ) - 1.0 );
}

double PowerLawUncorrelatedStep( const double mfp, const double a )
{
  return a * mfp * ( pow( RandomReal(), -1.0 / ( a ) ) - 1.0 );
}

double PowerLawCorrelatedTransmittance( const double s, const double a, const double l )
{
  return pow( (a*l)/(a*l + s), 1.0 + a);
}

double PowerLawUncorrelatedTransmittance( const double s, const double a, const double l )
{
  return pow( (a*l)/(a*l + s), a);
}

// random variate from Gamma distribution with paramter a > 1, from [Marsaglia and Tsang 2000]
double RandomGamma( const double a )
{
  double x,v,u;
  const double d = a - 1.0 / 3.0;
  const double c = 1.0 / sqrt( 9.0 * d );
  while( true )
  {
    do
    {
      x = RandomGauss();
      v = 1.0 + c * x;
    }
    while( v <= 0.0 );

    v = v * v * v;
    u = RandomReal();
    if( u < 1.0 - 0.0331 * (x*x)*(x*x) )
    {
      return d * v;
    }
    if( log(u) < 0.5 * x * x + d * ( 1.0 - v + log(v) ) )
    {
      return d * v;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////////
// Random Directions
//////////////////////////////////////////////////////////////////////////////////

Vector2 diskSample2D( const double radius )
{
  const double phi = RandomReal( 0.0, 2 * M_PI );
  return radius * sqrt( RandomReal() ) * Vector2( cos( phi ), sin( phi ) );
}

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

// sample a Lambertian direction about normal n:
Vector3 lambertDir(const Vector3& n)
{
  const Vector3 local(lambertDir());
  Vector3 x,y;
  buildOrthonormalBasis(x,y,n);
  return x * local.x + y * local.y + n * local.z;
}

Vector2 lambertDirFlatland()
{
    const double x = RandomReal();
    const double w = sqrt(2.0*x - x * x);
    const double w2 = RandomReal();
    return Vector2( sqrt( 1.0 - w * w ) * (w2 < 0.5 ? 1.0 : -1.0), w ); // y axis is the normal
}

Vector4 lambertDir4D()
{
  const double x = RandomReal(); 
  const double x4 = sqrt(1.0 - pow(1.0 - x, 0.6666666666666666 ) );
  return Vector4( sqrt( 1.0 - x4 * x4 ), x4, 0.0, 0.0 ); // y axis is the normal - watch out
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

// Rayleigh scattering
Vector3 RayleighDir( const Vector3& in )
{
    const double xi1 = RandomReal();
    const double w = (1 - Power(2 - 4*xi1 + Sqrt(5 + 16*(-1 + xi1)*xi1),0.6666666666666666))/
   Power(2 - 4*xi1 + Sqrt(5 + 16*(-1 + xi1)*xi1),0.3333333333333333);
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

// Lambertian Sphere Phase Function [Esposito and Lumme 1977]
Vector3 lambertSphereDir( const Vector3& in)
{
    const double xi1 = RandomReal();
    const double xi3 = RandomReal();
    const double xi4 = RandomReal();
    double w = -Sqrt(xi1*xi3) + Sqrt((-1 + xi1)*(-1 + xi3))*sin(2*Pi*xi4);
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