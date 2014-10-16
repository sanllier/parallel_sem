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
    checkres( MPI_Comm_size( MPI_COMM_WORLD, &procNum ) );
    checkres( MPI_Comm_rank( MPI_COMM_WORLD, &myID ) );
    
    //-----------------------------------------------------

    parparser arguments( argc, argv ); 
    const char* matAFile = arguments.get( "mata" ).value;
    const char* matBFile = arguments.get( "matb" ).value;
    const char* resFile = arguments.get( "res" ).value;
    if ( !matAFile || !matBFile || !resFile )
        throw "Invalid file name";

    const int xDim = arguments.get( "x" ).asInt(0);
    const int yDim = arguments.get( "y" ).asInt(0);
    const int zDim = arguments.get( "z" ).asInt(0);
    if ( xDim <= 0 || yDim <= 0 || zDim <= 0 )
        throw "Invalid processor cube dims";

    //----------------------------------------------------

    const int dims[3] = { xDim, yDim, zDim  };
    const int periods[3] = { 0, 0, 0 };
    MPI_Comm cubeComm;
    checkres( MPI_Cart_create( MPI_COMM_WORLD, 3, dims, periods, true, &cubeComm ) );

    int remainDims[3] = { 0, 0, 1 };
    MPI_Comm aDepthComm;
    checkres( MPI_Cart_sub( cubeComm, remainDims, &aDepthComm ) );

    remainDims[0] = 1; remainDims[1] = 0; remainDims[2] = 0;
    MPI_Comm bDepthComm;
    checkres( MPI_Cart_sub( cubeComm, remainDims, &bDepthComm ) );

    remainDims[0] = 0; remainDims[1] = 1;  remainDims[2] = 0;
    MPI_Comm cDepthComm;
    checkres( MPI_Cart_sub( cubeComm, remainDims, &cDepthComm ) );

    remainDims[0] = 1; remainDims[1] = 1; remainDims[2] = 0;
    MPI_Comm aSideComm;
    checkres( MPI_Cart_sub( cubeComm, remainDims, &aSideComm ) );

    remainDims[0] = 0; remainDims[1] = 1; remainDims[2] = 1;
    MPI_Comm bSideComm;
    checkres( MPI_Cart_sub( cubeComm, remainDims, &bSideComm ) );

    remainDims[0] = 1; remainDims[1] = 0; remainDims[2] = 1;
    MPI_Comm cSideComm;
    checkres( MPI_Cart_sub( cubeComm, remainDims, &cSideComm ) );

    int myADepthId = 0;
    checkres( MPI_Comm_rank( aDepthComm, &myADepthId ) );

    int myCDepthId = 0;
    checkres( MPI_Comm_rank( cDepthComm, &myCDepthId ) );

    int myCubeRank = 0;
    checkres( MPI_Comm_rank( cubeComm, &myCubeRank ) );

    int myASideRank = 0;
    checkres( MPI_Comm_rank( aSideComm, &myASideRank ) );

    int myBSideRank = 0;
    checkres( MPI_Comm_rank( bSideComm, &myBSideRank ) );

    int myCSideRank = 0;
    checkres( MPI_Comm_rank( cSideComm, &myCSideRank ) );

    int cubeCoords[3];
    checkres( MPI_Cart_coords( cubeComm, myCubeRank, 3, cubeCoords ) );

    int aSideCoords[2];
    checkres( MPI_Cart_coords( aSideComm, myASideRank, 2, aSideCoords ) );

    int bSideCoords[2];
    checkres( MPI_Cart_coords( bSideComm, myBSideRank, 2, bSideCoords ) );

    int cSideCoords[2];
    checkres( MPI_Cart_coords( cSideComm, myCSideRank, 2, cSideCoords ) );

    //----------------------------------------------------
  
    SHeader aHeader;
    SHeader bHeader;
    matrix< MATRIX_TYPE > aBlock;
    matrix< MATRIX_TYPE > bBlock;
    aBlock = parallelRead< MATRIX_TYPE, MPI_TYPE >( matAFile, myASideRank, xDim, yDim, aSideCoords, aSideComm, &aHeader );
    bBlock = parallelRead< MATRIX_TYPE, MPI_TYPE >( matBFile, myBSideRank, yDim, zDim, bSideCoords, bSideComm, &bHeader );

    if ( cubeCoords[2] != 0 )
        aBlock = matrix< MATRIX_TYPE >( aHeader.height / xDim, aHeader.width / yDim );
    if ( cubeCoords[0] != 0 )
        bBlock = matrix< MATRIX_TYPE >( bHeader.height / yDim, bHeader.width / zDim );

    MPI_Datatype ABLOCK = MPI_Type_vector_wrapper( aBlock.height(), aBlock.width(), aBlock.dataWidth(), MPI_TYPE );
    MPI_Datatype BBLOCK = MPI_Type_vector_wrapper( bBlock.height(), bBlock.width(), bBlock.dataWidth(), MPI_TYPE );
    checkres( MPI_Bcast( aBlock.raw(), 1, ABLOCK, ROOT_ID, aDepthComm ) );
    checkres( MPI_Bcast( bBlock.raw(), 1, BBLOCK, ROOT_ID, bDepthComm ) );

    //-----------------------------------------------------

    const double start = MPI_Wtime();

    matrix< MATRIX_TYPE >* resBlock = matrix_helper< MATRIX_TYPE >::mul( aBlock, bBlock );
    MPI_Datatype RESBLOCK = MPI_Type_vector_wrapper( resBlock->height(), resBlock->width(), resBlock->dataWidth(), MPI_TYPE );

    matrix< MATRIX_TYPE > tempBlock( resBlock->height(), resBlock->width() );
    if ( myCDepthId == ROOT_ID )
    {
        MPI_Status status;
        for ( int i = 1; i < yDim; ++i )
        {
            checkres( MPI_Recv( tempBlock.raw(), 1, RESBLOCK, MPI_ANY_SOURCE, MPI_ANY_TAG, cDepthComm, &status ) );
            resBlock->add( tempBlock );
        }
    }
    else
    {
        checkres( MPI_Send( resBlock->raw(), 1, RESBLOCK, ROOT_ID, 0, cDepthComm ) );
    }

    const double end = MPI_Wtime();
    MASTERPRINT( "TOTAL TIME: " << end - start << "\r\n" );

    parallelWrite< MATRIX_TYPE, MPI_TYPE >( resFile, myCSideRank, xDim, zDim, cSideCoords, cSideComm, *resBlock );

    delete resBlock;
    resBlock = 0;
    
    //-----------------------------------------------------

    checkres( MPI_Finalize() );
    return 0;
}
