#include <iostream>
#include <ctime>

#include "matrix\matrix.h"
#include "matrix\matrix_helper.h"
#include "matrix\matrix_serialization.h"

#include "parparser\parparser.h"

//----------------------------------------------------------

using namespace Matrix;

//----------------------------------------------------------

int main( int argc, char** argv )
{
    parparser arguments( argc, argv );
    const char* aFile = arguments.get( "af" ).value;
    const char* bFile = arguments.get( "bf" ).value;
    const char* resFile = arguments.get( "resf" ).value;

    clock_t start = clock();

    matrix<int> matA;
    matrix_serialization serializer;
    SHeader header;
    matA = *(matrix<int>*)serializer.readBinary( aFile, header );

    matrix<int> matB;
    matB = *(matrix<int>*)serializer.readBinary( bFile, header );

    matA.mul( matB );
    serializer.writeBinary( resFile, matA );

    clock_t end = clock();
    std::cout << "TOTAL TIME: " << end - start << "\r\n";

    return 0;
}
