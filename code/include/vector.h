#pragma once

//////////////////////////////////////////////////////////////////////////////////
// Vector2
//////////////////////////////////////////////////////////////////////////////////

struct Vector2{
    Vector2( const double in_x, const double in_y) 
    : x(in_x), y(in_y)
    {};

    Vector2() 
    : x(0.0), y(0.0)
    {};
    
  double x,y;

  Vector2 operator* (const double c ) const
  {
    return Vector2( c * this->x, c * this->y );
  };

  Vector2 operator/ (const double c ) const
  {
    return Vector2( this->x / c, this->y / c );
  };

  Vector2 operator+ ( Vector2 b ) const
  {
    return Vector2( b.x + this->x, b.y + this->y );
  };

  Vector2 operator- ( ) const
  {
    return Vector2( -this->x, -this->y );
  };

  Vector2 operator+= ( Vector2 b )
  {
    this->x += b.x;
    this->y += b.y;
    return *this;
  };

  Vector2 operator-= ( Vector2 b )
  {
    this->x -= b.x;
    this->y -= b.y;
    return *this;
  };
};

Vector2 operator*( const double c, const Vector2& v )
{
  return Vector2( c * v.x, c * v.y );
}

Vector2 operator-( const Vector2& a, const Vector2& b )
{
  return Vector2( a.x - b.x, a.y - b.y );
}

bool operator==( const Vector2& a, const Vector2& b )
{
  return a.x == b.x && a.y == b.y;
}

bool operator!=( const Vector2& a, const Vector2& b )
{
  return a.x != b.x || a.y != b.y;
}

double Dot( const Vector2& a, const Vector2& b )
{
  return a.x * b.x + a.y * b.y;
}

double Norm( const Vector2& v )
{
    return sqrt( v.x * v.x + v.y * v.y );
}

Vector2 Normalize( const Vector2& v )
{
  return v * ( 1.0 / Norm( v ) );
}

//////////////////////////////////////////////////////////////////////////////////
// Vector3
//////////////////////////////////////////////////////////////////////////////////

struct Vector3{
    Vector3( const double in_x, const double in_y, const double in_z) 
    : x(in_x), y(in_y), z(in_z)
    {};

    Vector3() 
    : x(0.0), y(0.0), z(0.0)
    {};
  double x,y,z;

  Vector3( const double * ptr )
  {
    x = ptr[0];
    y = ptr[1];
    z = ptr[2];
  }

  Vector3 operator* (const double c ) const
  {
    return Vector3( c * this->x, c * this->y, c * this->z );
  };

  Vector3 operator/ (const double c ) const
  {
    return Vector3( this->x / c, this->y / c, this->z / c );
  };

  Vector3 operator+ ( Vector3 b ) const
  {
    return Vector3( b.x + this->x, b.y + this->y, b.z + this->z );
  };

  Vector3 operator- ( ) const
  {
    return Vector3( -this->x, -this->y, -this->z );
  };

  Vector3 operator+= ( Vector3 b )
  {
    this->x += b.x;
    this->y += b.y;
    this->z += b.z;

    return *this;
  };

  Vector3 operator-= ( Vector3 b )
  {
    this->x -= b.x;
    this->y -= b.y;
    this->z -= b.z;

    return *this;
  };
};

Vector3 operator*( const double c, const Vector3& v )
{
  return Vector3( c * v.x, c * v.y, c * v.z );
}

Vector3 operator-( const Vector3& a, const Vector3& b )
{
  return Vector3( a.x - b.x, a.y - b.y, a.z - b.z );
}

bool operator==( const Vector3& a, const Vector3& b )
{
  return a.x == b.x && a.y == b.y && a.z == b.z;
}

bool operator!=( const Vector3& a, const Vector3& b )
{
  return a.x != b.x || a.y != b.y || a.z != b.z;
}

double Dot( const Vector3& a, const Vector3& b )
{
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

double Norm( const Vector3& v )
{
    return sqrt( v.x * v.x + v.y * v.y + v.z * v.z );
}

Vector3 Normalize( const Vector3& v )
{
  return v * (1.0 / Norm( v ) );
}

Vector3 Cross( const Vector3& a, const Vector3& b )
{
  return Vector3( -b.y * a.z + a.y * b.z, b.x * a.z - a.x * b.z, -b.x * a.y + a.x * b.y );
}

Vector3 reflect( const Vector3& in, const Vector3& n )
{
  return -in + 2.0 * Dot(in,n) * n;
}

//////////////////////////////////////////////////////////////////////////////////
// Vector4
//////////////////////////////////////////////////////////////////////////////////

struct Vector4{
    Vector4( const double in_x, const double in_y, const double in_z, const double in_w) 
    : x(in_x), y(in_y), z(in_z), w(in_w)
    {};
  double x,y,z,w;

  Vector4 operator* (const double c ) const
  {
    return Vector4( c * this->x, c * this->y, c * this->z, c * this->w );
  };

  Vector4 operator/ (const double c ) const
  {
    return Vector4( this->x / c, this->y / c, this->z / c, this->w / c );
  };

  Vector4 operator+ ( Vector4 b ) const
  {
    return Vector4( b.x + this->x, b.y + this->y, b.z + this->z, b.w + this->w );
  };

  Vector4 operator- ( ) const
  {
    return Vector4( -this->x, -this->y, -this->z, -this->w );
  };

  Vector4 operator+= ( Vector4 b )
  {
    this->x += b.x;
    this->y += b.y;
    this->z += b.z;
    this->w += b.w;

    return *this;
  };

  Vector4 operator-= ( Vector4 b )
  {
    this->x -= b.x;
    this->y -= b.y;
    this->z -= b.z;
    this->w -= b.w;

    return *this;
  };
};

Vector4 operator*( const double c, const Vector4& v )
{
  return Vector4( c * v.x, c * v.y, c * v.z, c * v.w );
}

Vector4 operator-( const Vector4& a, const Vector4& b )
{
  return Vector4( a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w );
}

bool operator==( const Vector4& a, const Vector4& b )
{
  return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}

bool operator!=( const Vector4& a, const Vector4& b )
{
  return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}

double Dot( const Vector4& a, const Vector4& b )
{
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

double Norm( const Vector4& v )
{
    return sqrt( v.x * v.x + v.y * v.y + v.z * v.z + v.w * v.w );
}

Vector4 Normalize( const Vector4& v )
{
  return v * (1.0 / Norm( v ) );
}

 // build orthonormal basis (Building an Orthonormal Basis from a 3D Unit Vector Without Normalization, [Frisvad2012])
 void buildOrthonormalBasis(Vector3& omega_1, Vector3& omega_2, const Vector3& omega_3)
{
	if(omega_3.z < -0.9999999f) 
	{
	   omega_1 = Vector3 ( 0.0f , -1.0f , 0.0f );
	   omega_2 = Vector3 ( -1.0f , 0.0f , 0.0f );
	} else {
	   const float a = 1.0f /(1.0f + omega_3.z );
	   const float b = -omega_3.x*omega_3 .y*a ;
	   omega_1 = Vector3 (1.0f - omega_3.x*omega_3. x*a , b , -omega_3.x );
	   omega_2 = Vector3 (b , 1.0f - omega_3.y*omega_3.y*a , -omega_3.y );
	}
}