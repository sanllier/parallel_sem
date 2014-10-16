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

    int myCDepthId = 0;
    checkres( MPI_Comm_rank( cDepthComm, &myCDepthId ) );

    int myCubeRank = 0;
    checkres( MPI_Comm_rank( cubeComm, &myCubeRank ) );

    int cubeCoords[3];
    checkres( MPI_Cart_coords( cubeComm, myCubeRank, 3, cubeCoords ) );

    //----------------------------------------------------
   
    const double start = MPI_Wtime();

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

        long matrixesDims[2] = { aBlock.height(), aBlock.width() };
        MPI_Bcast( matrixesDims, 2, MPI_LONG, ROOT_ID, MPI_COMM_WORLD );
        ABLOCK = MPI_Type_vector_wrapper( aBlock.height(), aBlock.width(), aBlock.dataWidth(), MPI_TYPE );
        MPI_Bcast( aBlock.raw(), 1, ABLOCK, ROOT_ID, aDepthComm );
    }
    else
    {
        long matrixesDims[2];
        MPI_Bcast( matrixesDims, 2, MPI_LONG, ROOT_ID, MPI_COMM_WORLD );

        ABLOCK = MPI_Type_vector_wrapper( matrixesDims[0], matrixesDims[1], matrixesDims[1], MPI_TYPE );

        aBlock = matrix< MATRIX_TYPE >( matrixesDims[0], matrixesDims[1] );
        MPI_Bcast( aBlock.raw(), 1, ABLOCK, ROOT_ID, aDepthComm );
    }
    
    if ( cubeCoords[0] == 0 )
    {
        matrix< MATRIX_TYPE > tempLine;

        SHeader bSrcMatHeader;
        SHeader bLineHeader;
        tempLine = deserializeLine< MATRIX_TYPE >( matBFile, cubeCoords[1], zDim, &bSrcMatHeader, &bLineHeader );
        long bBlockH = bSrcMatHeader.height / zDim;
        long bBlockW = bSrcMatHeader.width / yDim;
        bBlock.strongSubmatrix( tempLine, 0, bBlockW * cubeCoords[2], bBlockH, bBlockW );

        long matrixesDims[2] = { bBlock.height(), bBlock.width() };
        MPI_Bcast( matrixesDims, 2, MPI_LONG, ROOT_ID, MPI_COMM_WORLD );
        BBLOCK = MPI_Type_vector_wrapper( bBlock.height(), bBlock.width(), bBlock.dataWidth(), MPI_TYPE );
        MPI_Bcast( bBlock.raw(), 1, BBLOCK, ROOT_ID, bDepthComm );
    }
    else
    {
        long matrixesDims[2];
        MPI_Bcast( matrixesDims, 2, MPI_LONG, ROOT_ID, MPI_COMM_WORLD );

        BBLOCK = MPI_Type_vector_wrapper( matrixesDims[0], matrixesDims[1], matrixesDims[1], MPI_TYPE );

        bBlock = matrix< MATRIX_TYPE >( matrixesDims[0], matrixesDims[1] );
        MPI_Bcast( bBlock.raw(), 1, BBLOCK, ROOT_ID, bDepthComm );
    }

    //-----------------------------------------------------

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

    //-----------------------------------------------------

    if ( cubeCoords[1] == 0 )
        matrix_helper< MATRIX_TYPE >::print( *resBlock, std::cout );

    if ( cubeCoords[1] == 0 && cubeCoords[2] == 0 )
    {
        matrix< MATRIX_TYPE > resLine( resBlock->height(), resBlock->width() * zDim );
        resLine.insertSubmatrix( *resBlock, 0, 0 );

        MPI_Status status;
        for ( int i = 1; i < zDim; ++i )
        {
            checkres( MPI_Recv( resBlock->raw(), 1, RESBLOCK, MPI_ANY_SOURCE, MPI_ANY_TAG, cubeComm, &status ) );
            int coords[3];
            MPI_Cart_coords( cubeComm, status.MPI_SOURCE, 3, coords );
            resLine.insertSubmatrix( *resBlock, 0, coords[1] * resBlock->width() );
        }

        SHeader fullResHeader = { resLine.height() * xDim, resLine.width(), resLine.dataType(), resLine.matrixType() };
        SHeader lineResHeader = { resLine.height(), resLine.width(), resLine.dataType(), resLine.matrixType() };

        if ( cubeCoords[0] == 0 && cubeCoords[1] == 0 && cubeCoords[2] == 0 )
        {
            std::ofstream res( resFile, std::fstream::out | std::ios::binary );
            if ( !res.good() )
                return false;

            res.write( (const char*)&fullResHeader, sizeof( Matrix::SHeader ) );
            const char* raw = (char*)resLine.raw();
            res.write( raw, ELEMENT_SIZE_BY_TYPE[ lineResHeader.dataType ] * resLine.height() * resLine.width() );
            res.close();
        }

        MPI_Barrier( MPI_COMM_WORLD ); // CRAP

        for ( int i = 0; i < xDim; ++i )
        {
            if ( cubeCoords[0] == i && cubeCoords[1] == 0 && cubeCoords[2] == 0 )
            {
                std::fstream res( resFile, std::fstream::out | std::fstream::app | std::ios::binary );
                if ( !res.good() )
                    return false;

                const char* raw = (char*)resLine.raw();
                res.write( raw, ELEMENT_SIZE_BY_TYPE[ lineResHeader.dataType ] * resLine.height() * resLine.width() );
                res.close();
            }
            MPI_Barrier( MPI_COMM_WORLD ); // CRAP
        }
    }
    else if ( cubeCoords[1] == 0 && cubeCoords[2] > 0 )
    {
        int coords[3] = { cubeCoords[0], 0, 0 };
        int target;
        MPI_Comm_rank( cubeComm, &target );
        checkres( MPI_Send( resBlock->raw(), 1, RESBLOCK, target, 0, cubeComm ) );
        MPI_Barrier( MPI_COMM_WORLD ); // CRAP
        for ( int i = 0; i < xDim; ++i )
            MPI_Barrier( MPI_COMM_WORLD ); // CRAP
    }
    else
    {
        MPI_Barrier( MPI_COMM_WORLD ); // CRAP
        for ( int i = 0; i < xDim; ++i )
            MPI_Barrier( MPI_COMM_WORLD ); // CRAP
    }

    const double end = MPI_Wtime();
    MASTERPRINT( "TOTAL TIME: " << end - start << "\r\n" );

    delete resBlock;
    resBlock = 0;

    //-----------------------------------------------------

    checkres( MPI_Finalize() );
    return 0;
}
