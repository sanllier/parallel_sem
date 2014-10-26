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

    int dims[2] = { sqrtProcNum, sqrtProcNum  };
    int periods[2] = { 1, 1 };
    MPI_Comm gridComm;
    MPI_Cart_create( MPI_COMM_WORLD, 2, dims, periods, true, &gridComm );

    int coords[2];
    int myGridID = 0;
    MPI_Comm_rank( gridComm, &myGridID );
    MPI_Cart_coords( gridComm, myGridID, 2, coords );

    //----------------------------------------------------
    
    const double start = MPI_Wtime();

    SHeader aHeader;
    SHeader bHeader;
    matrix< MATRIX_TYPE > aBlock;
    matrix< MATRIX_TYPE > bBlock;
    aBlock = parallelRead< MATRIX_TYPE, MPI_TYPE >( matAFile, myID, sqrtProcNum, sqrtProcNum, coords, MPI_COMM_WORLD, &aHeader );
    bBlock = parallelRead< MATRIX_TYPE, MPI_TYPE >( matBFile, myID, sqrtProcNum, sqrtProcNum, coords, MPI_COMM_WORLD, &bHeader );

    if ( aHeader.width != bHeader.height )
        throw "Invalid matrices dimensions( non equal )";

    //-----------------------------------------------------

    int leftRank = 0;
    int rightRank = 0;
    int downRank = 0;
    int upRank = 0;
    checkres( MPI_Cart_shift( gridComm, 1, -1, &rightRank, &leftRank ) );
    checkres( MPI_Cart_shift( gridComm, 0, -1, &downRank, &upRank ) );

    MPI_Datatype SR_ABLOCK = MPI_Type_vector_wrapper( aBlock.height(), aBlock.width(), aBlock.dataWidth(), MPI_TYPE );
    MPI_Datatype SR_BBLOCK = MPI_Type_vector_wrapper( bBlock.height(), bBlock.width(), bBlock.dataWidth(), MPI_TYPE );

    int shiftSrc = 0;
    int shiftDst = 0;
    MPI_Status status;

    checkres( MPI_Cart_shift( gridComm, 1, -coords[0], &shiftSrc, &shiftDst ) );
    checkres( MPI_Sendrecv_replace( aBlock.shiftedRaw(), 1, SR_ABLOCK, shiftDst, SR_TAG, shiftSrc, SR_TAG, gridComm, &status ) );

    checkres( MPI_Cart_shift( gridComm, 0, -coords[1], &shiftSrc, &shiftDst ) );
    checkres( MPI_Sendrecv_replace( bBlock.shiftedRaw(), 1, SR_BBLOCK, shiftDst, SR_TAG, shiftSrc, SR_TAG, gridComm, &status ) );

    matrix< MATRIX_TYPE > resBlock( aBlock.height(), bBlock.width() );
    for ( int i = 0; i < dims[0]; ++i )
    {
        matrix_helper< MATRIX_TYPE >::mulsum( resBlock, aBlock, bBlock );
        checkres( MPI_Sendrecv_replace( aBlock.raw(), 1, SR_ABLOCK, leftRank, SR_TAG, rightRank, SR_TAG, gridComm, &status ) );
        checkres( MPI_Sendrecv_replace( bBlock.raw(), 1, SR_BBLOCK, upRank, SR_TAG, downRank, SR_TAG, gridComm, &status ) );
    }

    //-----------------------------------------------------

    parallelWrite< MATRIX_TYPE, MPI_TYPE >( resFile, myID, sqrtProcNum, sqrtProcNum, coords, MPI_COMM_WORLD, resBlock );

    const double end = MPI_Wtime();
    MASTERPRINT( "TOTAL TIME: " << end - start << "\r\n" );

    //-----------------------------------------------------

    checkres( MPI_Finalize() );
    return 0;
}
