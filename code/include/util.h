#pragma once

#include <cstdlib>
#include <sstream>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <cassert>

#define Sqrt sqrt
#define Power pow
#define Log log

//////////////////////////////////////////////////////////////////////////////////
// utilities
//////////////////////////////////////////////////////////////////////////////////

template <typename T>
T StringToNumber ( const std::string &Text )
{
    std::stringstream ss(Text);
    T result;
    return ss >> result ? result : 0;
}

// clamp a value into a specified range
double Clamp( const double x, const double min, const double max )
{
  return std::max( std::min( max, x ), min );
}

// map a continuous range [min,max] into dx sized intervals and return the integer index
// of an element x
size_t discreteMap( const double min, const double max, const double dx, const double x )
{
  const double clampx = Clamp( x, min, max );
  return size_t( floor( ( clampx - min ) / dx ) );
}