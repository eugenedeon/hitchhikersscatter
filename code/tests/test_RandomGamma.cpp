#include <util.h> // for discrete map
#include <random.h>

int main()
{
	for( int i = 0; i < 100000; ++i )
	{
		std::cout << RandomGamma(1.4) << " ";
	}
	std::cout << std::endl;
}