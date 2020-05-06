#include <util.h> // for discrete map
#include <random.h>

int main()
{
    for( size_t i = 0; i < 100000; i++ )
    {
    	std::cout << Chi3CorrelatedStep() << " ";
    }
    std::cout << std::endl;

	return 1;
}