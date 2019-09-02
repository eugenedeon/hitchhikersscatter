#include <fresnel.h>

int main()
{
    std::cout.precision(20);
	std::cout << "n = 1.0, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 ) << std::endl;
	std::cout << "n = 1.000000000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.000000000001 ) << std::endl;
	std::cout << "n = 1.00000000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.00000000001 ) << std::endl;
	std::cout << "n = 1.0000000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0000000001 ) << std::endl;
	std::cout << "n = 1.000000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.000000001 ) << std::endl;
	std::cout << "n = 1.00000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.00000001 ) << std::endl;
	std::cout << "n = 1.0000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0000001 ) << std::endl;
	std::cout << "n = 1.000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.000001 ) << std::endl;
	std::cout << "n = 1.00001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.00001 ) << std::endl;
	std::cout << "n = 1.0001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0001 ) << std::endl;
	std::cout << "n = 1.001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.001 ) << std::endl;
	std::cout << "n = 1.01, FR = " << SmoothDielectricHemisphericalAlbedo( 1.01 ) << std::endl;
	std::cout << "n = 1.1, FR = " << SmoothDielectricHemisphericalAlbedo( 1.1 ) << std::endl;
	std::cout << "n = 1.4, FR = " << SmoothDielectricHemisphericalAlbedo( 1.4 ) << std::endl;
	std::cout << "n = 2, FR = " << SmoothDielectricHemisphericalAlbedo( 2.0 ) << std::endl;
	std::cout << "n = 3, FR = " << SmoothDielectricHemisphericalAlbedo( 3.0 ) << std::endl;
	std::cout << "n = 10, FR = " << SmoothDielectricHemisphericalAlbedo( 10.0 ) << std::endl;
	std::cout << "n = 100, FR = " << SmoothDielectricHemisphericalAlbedo( 100.0 ) << std::endl;
	std::cout << "n = 1000, FR = " << SmoothDielectricHemisphericalAlbedo( 1000.0 ) << std::endl;
	std::cout << "n = 10000, FR = " << SmoothDielectricHemisphericalAlbedo( 10000.0 ) << std::endl;
	std::cout << "n = 100000, FR = " << SmoothDielectricHemisphericalAlbedo( 100000.0 ) << std::endl;

	std::cout << "n = 1/1.000000000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1.000000000001 ) << std::endl;
	std::cout << "n = 1/1.00000000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1.00000000001 ) << std::endl;
	std::cout << "n = 1/1.0000000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1.0000000001 ) << std::endl;
	std::cout << "n = 1/1.000000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1.000000001 ) << std::endl;
	std::cout << "n = 1/1.00000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1.00000001 ) << std::endl;
	std::cout << "n = 1/1.0000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1.0000001 ) << std::endl;
	std::cout << "n = 1/1.000001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1.000001 ) << std::endl;
	std::cout << "n = 1/1.00001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1.00001 ) << std::endl;
	std::cout << "n = 1/1.0001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1.0001 ) << std::endl;
	std::cout << "n = 1/1.001, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1.001 ) << std::endl;
	std::cout << "n = 1/1.01, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1.01 ) << std::endl;
	std::cout << "n = 1/1.1, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1.1 ) << std::endl;
	std::cout << "n = 1/1.4, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1.4 ) << std::endl;
	std::cout << "n = 1/2, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 2.0 ) << std::endl;
	std::cout << "n = 1/3, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 3.0 ) << std::endl;
	std::cout << "n = 1/10, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 10.0 ) << std::endl;
	std::cout << "n = 1/100, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 100.0 ) << std::endl;
	std::cout << "n = 1/1000, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 1000.0 ) << std::endl;
	std::cout << "n = 1/10000, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 10000.0 ) << std::endl;
	std::cout << "n = 1/100000, FR = " << SmoothDielectricHemisphericalAlbedo( 1.0 / 100000.0 ) << std::endl;
}