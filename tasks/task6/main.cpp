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
    // HOLY CRAP
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

    remainDims[0] = 1;  remainDims[1] = 1; remainDims[2] = 0;
    MPI_Comm aSideComm;
    checkres( MPI_Cart_sub( cubeComm, remainDims, &aSideComm ) );

    remainDims[0] = 0;  remainDims[1] = 1; remainDims[2] = 1;
    MPI_Comm bSideComm;
    checkres( MPI_Cart_sub( cubeComm, remainDims, &bSideComm ) );

    remainDims[0] = 1; remainDims[1] = 0; remainDims[2] = 1;
    MPI_Comm cSideComm;
    checkres( MPI_Cart_sub( cubeComm, remainDims, &cSideComm ) );

    int myCubeRank = 0;
    checkres( MPI_Comm_rank( cubeComm, &myCubeRank ) );

    int myADepthId = 0;
    checkres( MPI_Comm_rank( aDepthComm, &myADepthId ) );

    int myBDepthId = 0;
    checkres( MPI_Comm_rank( bDepthComm, &myBDepthId ) );

    int myCDepthId = 0;
    checkres( MPI_Comm_rank( cDepthComm, &myCDepthId ) );

    int myASideId = 0;
    checkres( MPI_Comm_rank( aSideComm, &myASideId ) );

    int myBSideId = 0;
    checkres( MPI_Comm_rank( bSideComm, &myBSideId ) );

    int cubeCoords[3];
    checkres( MPI_Cart_coords( cubeComm, myCubeRank, 3, cubeCoords ) );

    //----------------------------------------------------
   
    matrix< MATRIX_TYPE > aBlock;
    matrix< MATRIX_TYPE > bBlock;

    MPI_Datatype ABLOCK;
    MPI_Datatype BBLOCK;

    if ( cubeCoords[2] == 0 )
    {
        matrix< MATRIX_TYPE > tempLine;

        SHeader aSrcMatHeader;
        SHeader aLineHeader;
        tempLine = deserializeLine< MATRIX_TYPE >( matAFile, cubeCoords[0], xDim, &aSrcMatHeader, &aLineHeader );
        long aBlockH = aSrcMatHeader.height / xDim;
        long aBlockW = aSrcMatHeader.width / yDim;
        aBlock.strongSubmatrix( tempLine, 0, aBlockW * cubeCoords[1], aBlockH, aBlockW );

        SHeader bSrcMatHeader;
        SHeader bLineHeader;
        tempLine = deserializeLine< MATRIX_TYPE >( matBFile, cubeCoords[2], zDim, &bSrcMatHeader, &bLineHeader );
        long bBlockH = bSrcMatHeader.height / zDim;
        long bBlockW = bSrcMatHeader.width / yDim;
        bBlock.strongSubmatrix( tempLine, 0, bBlockW * cubeCoords[1], bBlockH, bBlockW );

        if ( aSrcMatHeader.width != bSrcMatHeader.height )
            throw "Invalid matrices dimensions( non equal )";

        long matrixesDims[4] = { aBlock.height(), aBlock.width(), bBlock.height(), bBlock.width() };
        MPI_Bcast( matrixesDims, 4, MPI_LONG, ROOT_ID, MPI_COMM_WORLD );

        ABLOCK = MPI_Type_vector_wrapper( aBlock.height(), aBlock.width(), aBlock.dataWidth(), MPI_TYPE );
        BBLOCK = MPI_Type_vector_wrapper( bBlock.height(), bBlock.width(), bBlock.dataWidth(), MPI_TYPE );

        MPI_Bcast( aBlock.raw(), 1, ABLOCK, ROOT_ID, aDepthComm );
        MPI_Bcast( bBlock.raw(), 1, BBLOCK, ROOT_ID, bDepthComm );
    }
    else
    {
        long matrixesDims[4];
        MPI_Bcast( matrixesDims, 4, MPI_LONG, ROOT_ID, MPI_COMM_WORLD );

        ABLOCK = MPI_Type_vector_wrapper( matrixesDims[0], matrixesDims[1], matrixesDims[1], MPI_TYPE );
        BBLOCK = MPI_Type_vector_wrapper( matrixesDims[2], matrixesDims[3], matrixesDims[3], MPI_TYPE );

        aBlock = matrix< MATRIX_TYPE >( matrixesDims[0], matrixesDims[1] );
        bBlock = matrix< MATRIX_TYPE >( matrixesDims[2], matrixesDims[3] );
        MPI_Bcast( aBlock.raw(), 1, ABLOCK, ROOT_ID, aDepthComm );
        MPI_Bcast( bBlock.raw(), 1, BBLOCK, ROOT_ID, bDepthComm );
    }
   
    //-----------------------------------------------------

    const double start = MPI_Wtime();

    matrix< MATRIX_TYPE >* resBlock = matrix_helper< MATRIX_TYPE >::mul( aBlock, bBlock );
    //MPI_Datatype RESBLOCK = MPI_Type_vector_wrapper( resBlock->height(), resBlock->width(), resBlock->dataWidth(), MPI_TYPE );

    //matrix< MATRIX_TYPE > tempBlock( resBlock->height(), resBlock->width() );
    //if ( myCDepthId == 0 )
    //{
    //    MPI_Status status;
    //    for ( int i = 0; i < yDim; ++i )
    //    {
    //        checkres( MPI_Recv( tempBlock.raw(), 1, RESBLOCK, MPI_ANY_SOURCE, MPI_ANY_TAG, cDepthComm, &status ) );
    //        resBlock->add( tempBlock );
    //    }
    //}

    matrix_helper< MATRIX_TYPE >::print( *resBlock, std::cout );

    delete resBlock;
    resBlock = 0;

    const double end = MPI_Wtime();
    MASTERPRINT( "TOTAL TIME: " << end - start << "\r\n" );

    //-----------------------------------------------------

    checkres( MPI_Finalize() );
    return 0;
}
