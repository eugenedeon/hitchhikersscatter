#include <Hfunction.h>

double psi3D( const double u, const double c )
{
    return c * 0.5;
}

double psi6D( const double u, const double c )
{
    return (8*c*Power(1 - Power(u,2),1.5))/(3.*Pi);
}

int main()
{
    std::cerr << Hfunction100( 1.0, 0.1, psi6D ) << std::endl;
}