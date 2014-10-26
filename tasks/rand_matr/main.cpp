#include <iostream>
#include <fstream>
#include <ctime>

#include "matrix/matrix.h"
#include "matrix/matrix_helper.h"
#include "matrix/matrix_serialization.h"

#include "parparser/parparser.h"

//----------------------------------------------------------

using namespace Matrix;

//---------------------------------------------------------

int main( int argc, char** argv )
{
    parparser arguments( argc, argv ); 
    const char* fileName = arguments.get( "f" ).value;
    if ( !fileName || !fileName[0] )
    {
        std::cout << "Invalid file name\r\n";
        return 0;
    }

    bool isPrint = arguments.get( "p" ).asBool( false );

    if ( !isPrint )
    {
        long height = arguments.get( "h" ).asLong( 4 );
        long width = arguments.get( "w" ).asLong( 4 );

        matrix<int> mat( height, width );
        srand( time(0) );
        matrix_helper<int>::fillRandom( mat );

        matrix_serialization serializer;
        serializer.writeBinary( fileName, mat );

        std::cout << height << " x " << width << " matrix was generated\r\n";  
    }
    else
    {
        matrix<int> mat;
        matrix_serialization serializer;
        SHeader header;
        mat = *(matrix<int>*)( serializer.readBinary( fileName, header ) );
        matrix_helper<int>::print( mat, std::cout );
    }

    return 0;
}
