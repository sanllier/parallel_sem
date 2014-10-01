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

enum
{
    SR_TAG
};

//----------------------------------------------------------

int main( int argc, char** argv )
{
    checkres( MPI_Init( &argc, &argv ) );

    int procNum = 0;
    int myID = 0;
    int sqrtProcNum = 0;
    checkres( MPI_Comm_size( MPI_COMM_WORLD, &procNum ) );
    checkres( MPI_Comm_rank( MPI_COMM_WORLD, &myID ) );
    sqrtProcNum = (int)( sqrt( procNum ) );  //CRAP
    
    //-----------------------------------------------------

    parparser arguments( argc, argv ); 
    const char* matAFile = arguments.get( "mata" ).value;
    const char* matBFile = arguments.get( "matb" ).value;
    const char* resFile = arguments.get( "res" ).value;
    if ( !matAFile || !matBFile || !resFile )
        throw "Invalid file name";

    //----------------------------------------------------

    const int dims[2] = { sqrtProcNum, sqrtProcNum  };
    const int periods[2] = { 1, 1 };
    MPI_Comm gridComm;
    MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, true, &gridComm );

    int coords[2];
    int myGridID = 0;
    MPI_Comm_rank( gridComm, &myGridID );
    MPI_Cart_coords( gridComm, myGridID, 2, coords );

    //----------------------------------------------------
    
    matrix< MATRIX_TYPE > aBlock;
    matrix< MATRIX_TYPE > bBlock;
    matrix< MATRIX_TYPE >* tempLine = new matrix< MATRIX_TYPE >;

    SHeader aSrcMatHeader;
    SHeader aLineHeader;
    *tempLine = deserializeLine< MATRIX_TYPE >( matAFile, coords[0], sqrtProcNum, &aSrcMatHeader, &aLineHeader );
    long aBlockH = aSrcMatHeader.height / sqrtProcNum;
    long aBlockW = aSrcMatHeader.width / sqrtProcNum;
    aBlock.strongSubmatrix( *tempLine, 0, aBlockW * coords[1], aBlockH, aBlockW );

    SHeader bSrcMatHeader;
    SHeader bLineHeader;
    *tempLine = deserializeLine< MATRIX_TYPE >( matBFile, coords[0], sqrtProcNum, &bSrcMatHeader, &bLineHeader );
    long bBlockH = bSrcMatHeader.height / sqrtProcNum;
    long bBlockW = bSrcMatHeader.width / sqrtProcNum;
    bBlock.strongSubmatrix( *tempLine, 0, bBlockW * coords[1], bBlockH, bBlockW );

    if ( aSrcMatHeader.width != bSrcMatHeader.height )
        throw "Invalid matrices dimensions( non equal )";

    delete tempLine;
    tempLine = 0;

    //-----------------------------------------------------

    int leftRank = 0;
    int rightRank = 0;
    int downRank = 0;
    int upRank = 0;
    MPI_Cart_shift( gridComm, 1, -1, &rightRank, &leftRank );
    MPI_Cart_shift( gridComm, 0, -1, &downRank, &upRank );

    MPI_Datatype SR_ABLOCK = MPI_Type_vector_wrapper( aBlock.height(), aBlock.width(), aBlock.dataWidth(), MPI_TYPE );
    MPI_Datatype SR_BBLOCK = MPI_Type_vector_wrapper( bBlock.height(), bBlock.width(), bBlock.dataWidth(), MPI_TYPE );

    int shiftSrc = 0;
    int shiftDst = 0;
    MPI_Status status;

    MPI_Cart_shift( gridComm, 1, -coords[0], &shiftSrc, &shiftDst );
    MPI_Sendrecv_replace( aBlock.shiftedRaw(), 1, SR_ABLOCK, shiftDst, SR_TAG, shiftSrc, SR_TAG, gridComm, &status );

    MPI_Cart_shift( gridComm, 0, -coords[1], &shiftSrc, &shiftDst );
    MPI_Sendrecv_replace( bBlock.shiftedRaw(), 1, SR_BBLOCK, shiftDst, SR_TAG, shiftSrc, SR_TAG, gridComm, &status );

    matrix< MATRIX_TYPE > resBlock( aBlock.height(), bBlock.width() );
    for ( int i = 0; i < dims[0]; ++i )
    {
        matrix_helper< MATRIX_TYPE >::mulsum( resBlock, aBlock, bBlock );
        MPI_Sendrecv_replace( aBlock.raw(), 1, SR_ABLOCK, leftRank, SR_TAG, rightRank, SR_TAG, gridComm, &status );
        MPI_Sendrecv_replace( bBlock.raw(), 1, SR_BBLOCK, upRank, SR_TAG, downRank, SR_TAG, gridComm, &status );
    }

    MPI_Cart_shift( gridComm, 1, coords[0], &shiftSrc, &shiftDst );
    MPI_Sendrecv_replace( aBlock.raw(), 1, SR_ABLOCK, shiftDst, SR_TAG, shiftSrc, SR_TAG, gridComm, &status );

    MPI_Cart_shift( gridComm, 0, coords[1], &shiftSrc, &shiftDst );
    MPI_Sendrecv_replace( bBlock.raw(), 1, SR_BBLOCK, shiftDst, SR_TAG, shiftSrc, SR_TAG, gridComm, &status );

    //-----------------------------------------------------

    MPI_Datatype SR_RESBLOCK = MPI_Type_vector_wrapper( resBlock.height(), resBlock.width(), resBlock.dataWidth(), MPI_TYPE );
    matrix< MATRIX_TYPE > resLine( resBlock.height(), resBlock.width() * sqrtProcNum );
    resLine.insertSubmatrix( resBlock, 0, myGridID * resBlock.width() );

    if ( coords[1] == 0 )
    {
        for ( int i = 1; i < sqrtProcNum; ++i ) 
        {
            MPI_Recv( resBlock.raw(), 1, SR_RESBLOCK, MPI_ANY_SOURCE, SR_TAG, gridComm, &status );
            int srcCoords[2];
            MPI_Cart_coords( gridComm, status.MPI_SOURCE, 2, coords );
            resLine.insertSubmatrix( resBlock, 0, srcCoords[1] * resBlock.width() );
        }
    }
    else
    {
        int dstRank = 0;
        int dstCoords[2] = { coords[0], 0 };
        MPI_Cart_rank( gridComm, dstCoords, &dstRank );
        MPI_Send( resBlock.raw(), 1, SR_RESBLOCK, dstRank, SR_TAG, gridComm );
    }

    //-----------------------------------------------------

    if ( coords[1] == 0 )
    {
        SHeader fullHeader = { resLine.height() * sqrtProcNum, resLine.width(), resLine.dataType(), resLine.matrixType() };
        SHeader lineHeader = { resLine.height(), resLine.width(), resLine.dataType(), resLine.matrixType() };
        const bool serializationRes = parallelMatrixSerialization( fullHeader, lineHeader, resLine, resFile );
        if ( !serializationRes )
            throw "Error while matrix serialization";
    }
    //-----------------------------------------------------

    checkres( MPI_Finalize() );
    return 0;
}
