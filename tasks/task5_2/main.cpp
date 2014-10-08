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

    const double start = MPI_Wtime();

    int leftRank = 0;
    int rightRank = 0;
    int downRank = 0;
    int upRank = 0;
    checkres( MPI_Cart_shift( gridComm, 1, -1, &rightRank, &leftRank ) );
    checkres( MPI_Cart_shift( gridComm, 0, -1, &downRank, &upRank ) );

    MPI_Datatype SR_ABLOCK = MPI_Type_vector_wrapper( aBlock.height(), aBlock.width(), aBlock.dataWidth(), MPI_TYPE );
    MPI_Datatype SR_BBLOCK = MPI_Type_vector_wrapper( bBlock.height(), bBlock.width(), bBlock.dataWidth(), MPI_TYPE );

    int remainDims[2] = { 0, 1 };
    MPI_Comm bcastComm;
    checkres( MPI_Cart_sub( gridComm, remainDims, &bcastComm ) );
        
    int shiftSrc = 0;
    int shiftDst = 0;
    MPI_Status status;

    matrix< MATRIX_TYPE > bcastedA( aBlock.height(), aBlock.width() );
    matrix< MATRIX_TYPE > workingBlock;
    matrix< MATRIX_TYPE > resBlock( aBlock.height(), bBlock.width() );
    for ( int i = 0; i < dims[0]; ++i )
    {
        if ( coords[1] == ( i + coords[0] ) % sqrtProcNum )
        {
            MPI_Bcast( aBlock.raw(), 1, SR_ABLOCK, coords[1], bcastComm );
            workingBlock = aBlock;
        }
        else
        {
            MPI_Bcast( bcastedA.raw(), 1, SR_ABLOCK, ( i + coords[0] ) % sqrtProcNum, bcastComm );
            workingBlock = bcastedA;
        }

        matrix_helper< MATRIX_TYPE >::mulsum( resBlock, workingBlock, bBlock );
        checkres( MPI_Sendrecv_replace( bBlock.raw(), 1, SR_BBLOCK, upRank, SR_TAG, downRank, SR_TAG, gridComm, &status ) );
    }
   
    //-----------------------------------------------------

    remainDims[0] = 1;
    remainDims[1] = 0;
    MPI_Comm outputComm;
    MPI_Cart_sub( gridComm, remainDims, &outputComm );

    MPI_Datatype SR_RESBLOCK = MPI_Type_vector_wrapper( resBlock.height(), resBlock.width(), resBlock.dataWidth(), MPI_TYPE );
    if ( coords[1] == 0 )
    {        
        matrix< MATRIX_TYPE > resLine( resBlock.height(), resBlock.width() * sqrtProcNum );
        resLine.insertSubmatrix( resBlock, 0, 0 );

        int srcCoords[2] = { 0, 0 };
        for ( int i = 1; i < sqrtProcNum; ++i ) 
        {
            checkres( MPI_Recv( resBlock.raw(), 1, SR_RESBLOCK, MPI_ANY_SOURCE, SR_TAG, gridComm, &status ) );
            checkres( MPI_Cart_coords( gridComm, status.MPI_SOURCE, 2, srcCoords ) );
            resLine.insertSubmatrix( resBlock, 0, srcCoords[1] * resBlock.width() );
        }

        // CRAP
        SHeader fullHeader = { resLine.height() * sqrtProcNum, resLine.width(), resLine.dataType(), resLine.matrixType() };
        SHeader lineHeader = { resLine.height(), resLine.width(), resLine.dataType(), resLine.matrixType() };
                  
        for ( int i = 0; i < procNum; ++i )
        {
            if ( i == myID )
            {
                if ( coords[0] == 0 )
                {
                    std::ofstream res( resFile, std::fstream::out | std::ios::binary );
                    if ( !res.good() )
                        return false;

                    res.write( (const char*)&fullHeader, sizeof( Matrix::SHeader ) );
                    const char* raw = (char*)resLine.raw();
                    res.write( raw, ELEMENT_SIZE_BY_TYPE[ lineHeader.dataType ] * resLine.height() * resLine.width() );
                    res.close();
                }
                else
                {
                    std::fstream res( resFile, std::fstream::out | std::fstream::app | std::ios::binary );
                    if ( !res.good() )
                        return false;

                    const char* raw = (char*)resLine.raw();
                    res.write( raw, ELEMENT_SIZE_BY_TYPE[ lineHeader.dataType ] * resLine.height() * resLine.width() );
                    res.close();
                }
            }
            MPI_Barrier( outputComm );
        }       
    }
    else
    {
        int dstRank = 0;
        int dstCoords[2] = { coords[0], 0 };
        checkres( MPI_Cart_rank( gridComm, dstCoords, &dstRank ) );
        checkres( MPI_Send( resBlock.raw(), 1, SR_RESBLOCK, dstRank, SR_TAG, gridComm ) );
    }

    const double end = MPI_Wtime();
    MASTERPRINT( "TOTAL TIME: " << end - start << "\r\n" );

    //-----------------------------------------------------

    checkres( MPI_Finalize() );
    return 0;
}
