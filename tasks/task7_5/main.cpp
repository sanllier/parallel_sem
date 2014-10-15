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
    
    SHeader aHeader;
    SHeader bHeader;
    if ( myID == ROOT_ID )
    {
        std::ifstream inpAFile( matAFile, std::ios::binary );
        std::ifstream inpBFile( matBFile, std::ios::binary );
        if ( !inpAFile.good() || !inpBFile.good() )
            throw "Error while read matrix file";

        inpAFile.read( (char *)( &aHeader ), sizeof( aHeader ) );
        inpBFile.read( (char *)( &bHeader ), sizeof( bHeader ) );

        inpAFile.close();
        inpBFile.close();

        MPI_Bcast( &aHeader, sizeof( aHeader ), MPI_CHAR, ROOT_ID, MPI_COMM_WORLD );
        MPI_Bcast( &bHeader, sizeof( bHeader ), MPI_CHAR, ROOT_ID, MPI_COMM_WORLD );
    }
    else
    {
        MPI_Bcast( &aHeader, sizeof( aHeader ), MPI_CHAR, ROOT_ID, MPI_COMM_WORLD );
        MPI_Bcast( &bHeader, sizeof( bHeader ), MPI_CHAR, ROOT_ID, MPI_COMM_WORLD );
    }

    if ( aHeader.width != bHeader.height )
        throw "Invalid matrices dimensions( non equal )";

    matrix< MATRIX_TYPE > aBlock( aHeader.height / sqrtProcNum, aHeader.width / sqrtProcNum );
    matrix< MATRIX_TYPE > bBlock( bHeader.height / sqrtProcNum, bHeader.width / sqrtProcNum );

    int gsizes[]   = { aHeader.height, aHeader.width };
    int distribs[] = { MPI_DISTRIBUTE_BLOCK, MPI_DISTRIBUTE_BLOCK };
    int dargs[]    = { MPI_DISTRIBUTE_DFLT_DARG, MPI_DISTRIBUTE_DFLT_DARG };
    int psizes[]   = { sqrtProcNum, sqrtProcNum };

    MPI_Datatype darrayType;
    MPI_File aMatFile;
    MPI_Status status;
    checkres( MPI_Type_create_darray( procNum, myID, 2, gsizes, distribs, dargs, psizes, MPI_ORDER_C, MPI_TYPE, &darrayType ) );
    checkres( MPI_Type_commit( &darrayType ) );
    checkres( MPI_File_open( MPI_COMM_WORLD, matAFile, MPI_MODE_RDONLY,  MPI_INFO_NULL, &aMatFile ) );
   // checkres( MPI_File_set_view( aMatFile, sizeof( aHeader ), MPI_TYPE, darrayType, "native", MPI_INFO_NULL) );
/*    MPI_File_read_all( aMatFile, aBlock.raw(), aBlock.height() * aBlock.width(), MPI_TYPE, &status );
    MPI_File_close( &aMatFile );

    matrix_helper< MATRIX_TYPE >::print( aBlock, std::cout );


/*



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

    checkres( MPI_Cart_shift( gridComm, 1, coords[0], &shiftSrc, &shiftDst ) );
    checkres( MPI_Sendrecv_replace( aBlock.raw(), 1, SR_ABLOCK, shiftDst, SR_TAG, shiftSrc, SR_TAG, gridComm, &status ) );

    checkres( MPI_Cart_shift( gridComm, 0, coords[1], &shiftSrc, &shiftDst ) );
    checkres( MPI_Sendrecv_replace( bBlock.raw(), 1, SR_BBLOCK, shiftDst, SR_TAG, shiftSrc, SR_TAG, gridComm, &status ) );

    //-----------------------------------------------------

    int remainDims[2] = { 1, 0 };
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
    */
    //-----------------------------------------------------

    checkres( MPI_Finalize() );
    return 0;
}
