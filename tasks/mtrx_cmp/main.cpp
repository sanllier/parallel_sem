#include <iostream>
#include <fstream>
#include "parparser\parparser.h"

static const size_t BUF_SIZE = 128;

int main( int argc, char** argv )
{
    parparser arguments( argc, argv );

    const char* aFile = arguments.get( "af" ).value;
    const char* bFile = arguments.get( "bf" ).value;
    if ( !aFile || !bFile )
    {
        std::cout << "Invalid file name\r\n";
        return 0;
    }

    std::ifstream aFStr( aFile, std::ios::binary );
    std::ifstream bFStr( bFile, std::ios::binary );
    if ( !aFStr.good() || !bFStr.good() )
    {
        aFStr.close();
        bFStr.close();
        std::cout << "Invalud file\r\n";
        return 0;
    }

    char aBuf[ BUF_SIZE ];
    char bBuf[ BUF_SIZE ];
    while ( !aFStr.eof() && !bFStr.eof() )
    {
        aFStr.read( aBuf, BUF_SIZE );
        bFStr.read( bBuf, BUF_SIZE );
        for ( size_t i = 0; i < BUF_SIZE; ++i )
        {
            if ( aBuf[i] != bBuf[i] )
            {
                std::cout << "Different\r\n";
                aFStr.close();
                bFStr.close();
                return 0;
            }
        }
    }

    if ( !aFStr.eof() || !bFStr.eof() )
    {
        std::cout << "Different by eof\r\n";
        aFStr.close();
        bFStr.close();
        return 0;
    }

    std::cout << "OK\r\n";
    aFStr.close();
    bFStr.close();
    return 0;
}
