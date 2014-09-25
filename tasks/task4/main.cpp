#include <iostream>
#include <fstream>
#include <vector>

#include "helper.h"

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
    void* aLine = deserializeLine( matFile, myID, procNum, &srcMatHeader );
    if ( srcMatHeader.height != srcMatHeader.width )
        throw "Incorrect matrix dimensions";

    double start = MPI_Wtime();

    int blockSize = srcMatHeader.height / procNum;
    matrix<int> result( aLine.height(), aLine.width() );
    MPI_Datatype SEND_MATRIX_BLOCK;
    MPI_Datatype RECV_MATRIX_BLOCK;
    MPI_Request sendRequest;
    std::vector< MPI_Request > recvRequests;
    std::vector< MPI_Status >  recvStatuses;

    for ( int _ = 1; _ < power; ++_  )
    {
        recvRequests.clear();
        recvStatuses.clear();

        for ( int iter = 0; iter < procNum; ++iter )
        {
            matrix<int> myBlock = getBlock( aLine, procNum, iter );

            if ( iter == 0 )
            {
                checkres( MPI_Type_vector( myBlock.height(), myBlock.width(), myBlock.dataWidth(), MPI_INT, &SEND_MATRIX_BLOCK ) );
                checkres( MPI_Type_commit( &SEND_MATRIX_BLOCK ) );
            }

            for ( int i = 0; i < procNum; ++i )
                if ( i != myID )
                    checkres( MPI_Isend( myBlock.shiftedRaw(), 1, SEND_MATRIX_BLOCK, i, BLOCK_BCAST, MPI_COMM_WORLD, &sendRequest) );

            matrix<int> bTranspLine( aLine.width(), blockWidth );
            for ( long i = 0; i < aLine.height(); ++i )
                for ( long q = 0; q < blockWidth; ++q )
                    bTranspLine.at( i + myID * aLine.height(), q ) = myBlock.at( i, q );
        
            if ( iter == 0 )
            {
                checkres( MPI_Type_vector( myBlock.height(), myBlock.width(), bTranspLine.dataWidth(), MPI_INT, &RECV_MATRIX_BLOCK ) );
                checkres( MPI_Type_commit( &RECV_MATRIX_BLOCK ) );
            }

            for ( int i = 0; i < procNum; ++i )
            {
                if ( i != myID )
                {
                    char* dataPos = (char *)bTranspLine.raw() + i * aLine.height() * blockWidth * ELEMENT_SIZE_BY_TYPE[ bTranspLine.dataType() ];
                    recvRequests.push_back( MPI_Request() );
                    recvStatuses.push_back( MPI_Status() );
                    checkres( MPI_Irecv( dataPos, 1, RECV_MATRIX_BLOCK, i, BLOCK_BCAST, MPI_COMM_WORLD, &recvRequests.back() ) );
                }
            }
            checkres( MPI_Waitall( recvRequests.size(), recvRequests.data(), recvStatuses.data() ) );

            if ( myID == 1 )
                matrix_helper<int>::print( bTranspLine, std::cout );
            if ( myID == 1 )
                matrix_helper<int>::print( aLine, std::cout );

            matrix<int>* resBlock = matrix_helper<int>::mul( aLine, bTranspLine );
            long targetCol = iter * blockWidth;
            for ( long i = 0; i < resBlock->height(); ++i )
                for ( long q = 0; q < resBlock->width(); ++q )
                    result.at( i, q + targetCol ) = resBlock->at( i, q );    

            if ( myID == 1 )
                matrix_helper<int>::print( result, std::cout );

            aLine = result;
            delete resBlock;
        }
    }

    double end = MPI_Wtime();
    MASTERPRINT( "TOTAL TIME: " << end - start << "\r\n" );

    //-----------------------------------------------------
    if ( power <= 1 )
        result = aLine;

    int id = 0;
    while ( id < procNum )
    {
        if ( myID == id )
        {
            if ( myID == ROOT_ID )
            {
                SHeader resHeader = { resH, resW, header.dataType, header.matrixType };
                std::ofstream res( resFile, std::ios::binary );
                res.write( (const char*)&resHeader, sizeof( SHeader ) );
                const char* raw = (char*)result.raw();
                res.write( raw, ELEMENT_SIZE_BY_TYPE[ header.dataType ] * result.height() * result.width() );
                res.close();
            }
            else
            {
                std::fstream res;
                res.open ( resFile, std::fstream::out | std::fstream::app | std::ios::binary );
                const char* raw = (char*)result.raw();
                res.write( raw, ELEMENT_SIZE_BY_TYPE[ header.dataType ] * result.height() * result.width() );
                res.close();
            }
        }
        ++id;
        MPI_Barrier( MPI_COMM_WORLD );
    }

    //-----------------------------------------------------
    
    checkres( MPI_Finalize() );
    return 0;
}
