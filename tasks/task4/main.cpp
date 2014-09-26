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
    int power = arguments.get("n").asInt(1);
    const char* matFile = arguments.get("f").value;
    const char* resFile = arguments.get("res").value;
    if ( !matFile || !matFile[0] || !resFile || !resFile[0] )
        throw "Invalid file name";

    //-----------------------------------------------------

    SHeader srcMatHeader;
    SHeader lineHeader;
    matrix< MATRIX_TYPE > aLine = deserializeLine< MATRIX_TYPE >( matFile, myID, procNum, &srcMatHeader, &lineHeader );
    if ( srcMatHeader.height != srcMatHeader.width )
        throw "Incorrect matrix dimensions";

    matrix< MATRIX_TYPE > aLineBackup( aLine );
       
    //-----------------------------------------------------

    int blockSize = srcMatHeader.width / procNum;
    matrix< MATRIX_TYPE > bTranspLine( aLine.width(), aLine.height() );
    matrix< MATRIX_TYPE > result( aLine.height(), aLine.width() );
    matrix< MATRIX_TYPE > resBlock( blockSize, blockSize );

    MPI_Datatype SEND_MATRIX_BLOCK = MPI_Type_vector_wrapper( blockSize, blockSize, aLine.dataWidth(), MPI_TYPE );
    MPI_Datatype RECV_MATRIX_BLOCK = MPI_Type_vector_wrapper( blockSize, blockSize, bTranspLine.dataWidth(), MPI_TYPE );

    double start = MPI_Wtime();

    for ( int _ = 1; _ < power; ++_  )
    {
        for ( int iter = 0; iter < procNum; ++iter )
        {
            for ( int i = 0; i < procNum; ++i )
            {
                if ( i == myID )
                {
                    matrix< MATRIX_TYPE > myBlock = getBlock( aLineBackup, procNum, iter ); 
                    MPI_Bcast( myBlock.shiftedRaw(), 1, SEND_MATRIX_BLOCK, i, MPI_COMM_WORLD );
                    bTranspLine.insertSubmatrix( myBlock, myID * blockSize, 0 );
                }
                else
                {            
                    char* dataPos = (char *)bTranspLine.raw() + i * blockSize * blockSize * ELEMENT_SIZE_BY_TYPE[ bTranspLine.dataType() ];
                    MPI_Bcast( dataPos, 1, RECV_MATRIX_BLOCK, i, MPI_COMM_WORLD );   
                }
            }

            matrix_helper< MATRIX_TYPE >::rmul( resBlock, aLine, bTranspLine );
            result.insertSubmatrix( resBlock, 0, iter * blockSize );
        }
        aLine = result;
    }

    double end = MPI_Wtime();
    MASTERPRINT( "TOTAL TIME: " << end - start << "\r\n" );

    //-----------------------------------------------------

    const bool serializationRes = parallelMatrixSerialization( srcMatHeader, lineHeader, result, resFile );
    if ( !serializationRes )
        throw "Error while matrix serialization";

    //-----------------------------------------------------
    
    checkres( MPI_Finalize() );
    return 0;
}
