#include "../GENHalfspaceFlatland.h"

int main( int argc, char** argv )
{
    if( argc != 5 )
    {
        std::cout << "usage: albedoEscapeProb c dz numz numsamples \n";
        exit( -1 );
    }

    double c = StringToNumber<double>( std::string( argv[1] ) );
    double dz = StringToNumber<double>( std::string( argv[2] ) );
    size_t numz = StringToNumber<double>( std::string( argv[3] ) );
    size_t numsamples = StringToNumber<size_t>( std::string( argv[4] ) );

    std::cout << "Albedo problem Flatland Escape Probability, c = " << c << std::endl;

    // simulate and estimate
    for( double minz = 0.0; minz < numz * dz; minz += dz )
    {
        const double maxz = minz + dz;
        size_t escapeCount = 0;

        for( size_t i = 0; i < numsamples; ++i )
        {
            const double startz = RandomReal( minz, maxz );

            Vector2 pos( 0.0, startz );
            Vector2 dir = isotropicDirFlatland();
            const double step = -log( RandomReal() );
            pos += dir * step;

            bool absorb = false;

            while( pos.y > 0.0 )
            {
                if( RandomReal() > c )
                {
                    absorb = true;
                    break;
                }

                Vector2 dir = isotropicDirFlatland();
                const double step = -log( RandomReal() );
                pos += dir * step;
            }

            if( !absorb )
                escapeCount++;
        }

        std::cout << double(escapeCount ) / double( numsamples ) << " ";
    }

    std::cout << std::endl;

    return 0;
}