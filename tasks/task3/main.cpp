#include <iostream>
#include <fstream>
#include <vector>

#include "helper.h"

//----------------------------------------------------------

// TODO: fix
#ifndef MATRIX_TYPE
    #define MATRIX_TYPE int
    #define MPI_TYPE MPI_INT 
#endif
//----------------------------------------------------------

using namespace Matrix;

//----------------------------------------------------------

int main( int argc, char** argv )
{
    checkres( MPI_Init( &argc, &argv ) );

    int procNum = 0;
    int myID = 0;
    checkres( MPI_Comm_size( MPI_COMM_WORLD, &procNum ) );
    checkres( MPI_Comm_rank( MPI_COMM_WORLD, &myID ) );

    //-----------------------------------------------------

    parparser arguments( argc, argv ); 
    const char* matAFile = arguments.get( "mata" ).value;
    const char* matBFile = arguments.get( "matb" ).value;
    const char* resFile = arguments.get( "res" ).value;
    if ( !matAFile || !matBFile || !resFile )
        throw "Invalid file name";

    //----------------------------------------------------
    
    SHeader aSrcMatHeader;
    SHeader aLineHeader;
    matrix< MATRIX_TYPE > aLine = deserializeLine< MATRIX_TYPE >( matAFile, myID, procNum, &aSrcMatHeader, &aLineHeader );

    SHeader bSrcMatHeader;
    SHeader bLineHeader;
    matrix< MATRIX_TYPE > bLine = deserializeLine< MATRIX_TYPE >( matBFile, myID, procNum, &bSrcMatHeader, &bLineHeader );

    if ( aSrcMatHeader.width != bSrcMatHeader.height )
        throw "Invalid matrices dimensions( non equal )";

    //----------------------------------------------------
 
    int blockWidth = bLine.width() / procNum;
    matrix< MATRIX_TYPE > bTranspLine( bSrcMatHeader.height, blockWidth );
    matrix< MATRIX_TYPE > resBlock( aLine.height(), bTranspLine.width() );
    matrix< MATRIX_TYPE > result( resBlock.height(), resBlock.width() * procNum );

    MPI_Datatype SEND_MATRIX_BLOCK = MPI_Type_vector_wrapper( bLine.height(), blockWidth, bLine.dataWidth(), MPI_TYPE );
    MPI_Datatype RECV_MATRIX_BLOCK = MPI_Type_vector_wrapper( bLine.height(), blockWidth, bTranspLine.dataWidth(), MPI_TYPE );

    double start = MPI_Wtime();

    for ( int iter = 0; iter < procNum; ++iter )
    {
        for ( int i = 0; i < procNum; ++i )
        {
            if ( i == myID )
            {
                matrix< MATRIX_TYPE > myBlock = getBlock( bLine, procNum, iter );
                MPI_Bcast( myBlock.shiftedRaw(), 1, SEND_MATRIX_BLOCK, i, MPI_COMM_WORLD );
                bTranspLine.insertSubmatrix( myBlock, myID * myBlock.height(), 0 );
            }
            else
            {            
                char* dataPos = (char *)bTranspLine.raw() + i * aLine.height() * blockWidth * ELEMENT_SIZE_BY_TYPE[ bTranspLine.dataType() ];
                MPI_Bcast( dataPos, 1, RECV_MATRIX_BLOCK, i, MPI_COMM_WORLD );   
            }
        }

        matrix_helper< MATRIX_TYPE >::rmul( resBlock, aLine, bTranspLine );
        result.insertSubmatrix( resBlock, 0, iter * blockWidth );
    }

    double end = MPI_Wtime();
    MASTERPRINT( "TOTAL TIME: " << end - start << "\r\n" );

    //-----------------------------------------------------

    SHeader fullResHeader = { aSrcMatHeader.height, bSrcMatHeader.width, aLine.dataType(), aLine.matrixType() };
    SHeader partResHeader = { result.height(), result.width(), aLine.dataType(), aLine.matrixType() };
    const bool serializationRes = parallelMatrixSerialization( fullResHeader, partResHeader, result, resFile );
    if ( !serializationRes )
        throw "Error while matrix serialization";

    //-----------------------------------------------------
    
    checkres( MPI_Finalize() );
    return 0;
}
