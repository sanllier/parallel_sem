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

void checkRoutine( parparser& arguments )
{
    int power = arguments.get("n").asInt(1);    
    const char* matFile = arguments.get("f").value;
    const char* resFile = arguments.get("res").value;
    if ( !matFile || !matFile[0] || !resFile || !resFile[0] )
        throw "Invalid file name";

    //----------------------------------------------------

    SHeader header;
    matrix_serialization serializer;
    matrix< MATRIX_TYPE >* workMatrix = ( matrix< MATRIX_TYPE >* )serializer.readBinary( matFile, header );
    matrix< MATRIX_TYPE > backUp( *workMatrix ); 

    for ( int i = 1; i < power; ++i )
        workMatrix->mul( backUp );

    serializer.writeBinary( resFile, *workMatrix );

    delete workMatrix;
    workMatrix = 0;

    //----------------------------------------------------

    exit(0);
}

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
    if ( arguments.get("check").asBool( false ) )
    {
        checkres( MPI_Finalize() );
        checkRoutine( arguments );
    }

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

    int blockSize = aLineBackup.width() / procNum;
    matrix< MATRIX_TYPE > bTranspLine( srcMatHeader.height, blockSize );
    matrix< MATRIX_TYPE > resBlock( aLine.height(), bTranspLine.width() );
    matrix< MATRIX_TYPE > result( resBlock.height(), resBlock.width() * procNum );
    matrix< MATRIX_TYPE >* workLine = &aLineBackup;

    MPI_Datatype SEND_MATRIX_BLOCK = MPI_Type_vector_wrapper( workLine->height(), blockSize, workLine->dataWidth(), MPI_TYPE );
    MPI_Datatype RECV_MATRIX_BLOCK = MPI_Type_vector_wrapper( workLine->height(), blockSize, bTranspLine.dataWidth(), MPI_TYPE );

    double start = MPI_Wtime();

    for ( int curPow = 1; curPow < power; ++curPow  )
    {
        for ( int iter = 0; iter < procNum; ++iter )
        {
            for ( int i = 0; i < procNum; ++i )
            {
                if ( i == myID )
                {
                    matrix< MATRIX_TYPE > myBlock = getBlock( *workLine, procNum, iter ); 
                    MPI_Bcast( myBlock.shiftedRaw(), 1, SEND_MATRIX_BLOCK, i, MPI_COMM_WORLD );
                    bTranspLine.insertSubmatrix( myBlock, myID * blockSize, 0 );
                }
                else
                {        
                    char* dataPos = (char *)bTranspLine.raw() + i * aLine.height() * blockSize * ELEMENT_SIZE_BY_TYPE[ bTranspLine.dataType() ];
                    MPI_Bcast( dataPos, 1, RECV_MATRIX_BLOCK, i, MPI_COMM_WORLD );   
                }
            }

            matrix_helper< MATRIX_TYPE >::rmul( resBlock, aLine, bTranspLine );
            result.insertSubmatrix( resBlock, 0, iter * blockSize );
        }
        aLine.equalStrongCopy( result );

        if ( power / ( curPow + 1 ) >= 2 )
        {
            workLine = &aLine;
            curPow = ( ( curPow + 1 ) << 1 ) - 2;
        }
        else
        {
            workLine = &aLineBackup;
        }
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
