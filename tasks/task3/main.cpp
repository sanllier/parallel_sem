#include <iostream>
#include <fstream>
#include <vector>

#include "matrix\matrix.h"
#include "matrix\matrix_helper.h"
#include "matrix\matrix_serialization.h"
//#include "matrix\matrix_mpi.h"
#include "mpi.h"
#include "parparser\parparser.h"

//----------------------------------------------------------

using namespace Matrix;

//----------------------------------------------------------

#define MASTERPRINT( _x_ ) if ( myID == ROOT_ID ) std::cout << _x_ ; std::cout.flush()

enum
{
    ROOT_ID = 0,
    BLOCK_BCAST
};

static const size_t BUF_SIZE = sizeof( SHeader );
static char buffer[ BUF_SIZE ];

inline void checkres( int res )
{
    if ( res != MPI_SUCCESS )
        throw res;
}

inline size_t fileSize( std::ifstream& iFstr )
{
    if ( !iFstr.good() )
        return 0;

    std::streampos beginPos = iFstr.tellg();
    iFstr.seekg( 0, std::ios::end );
    std::streampos endPos = iFstr.tellg();
    size_t size = (size_t)( endPos - beginPos );
    iFstr.seekg( beginPos );
    return size;
}

template< class T >
inline matrix<T> getBlock( matrix<T>& mat, int procNum, int block )
{
    long blockWidth = mat.width() / procNum;
    matrix<T> subMat( mat.height(), blockWidth );
    subMat.weakSubmatrix( mat, 0, block * blockWidth, mat.height(), blockWidth );
    return subMat;
}

//---------------------------------------------------------

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

    std::ifstream aFile( matAFile, std::ios::binary );
    std::ifstream bFile( matBFile, std::ios::binary );
    if ( !aFile.good() || !bFile.good() )
        throw "Invalid file";

    //----------------------------------------------------
    long resH = 0;
    long resW = 0;

    matrix_serialization serializer;
    size_t size = fileSize( aFile );
    aFile.read( buffer, BUF_SIZE );
    SHeader header = serializer.deserializeStart( BUF_SIZE, buffer );
    resH = header.height;
    if ( header.height % procNum || header.width % procNum )
        throw "Invalid matrix size";

    header.height /= procNum;
    int myShift = sizeof( SHeader ) + header.height * myID * header.width * ELEMENT_SIZE_BY_TYPE[ header.dataType ];
    aFile.seekg( myShift );
    size = header.height * header.width * ELEMENT_SIZE_BY_TYPE[ header.dataType ];
    matrix<int> aLine( header.height, header.width );

    while ( size > 0 )
    {
        size_t readSize = size < BUF_SIZE ? size : BUF_SIZE;
        aFile.read( buffer, readSize );
        serializer.deserializeStep( buffer, readSize, aLine );
        size -= readSize;        
    }
    serializer.deserializeStop();

    size = fileSize( bFile );
    bFile.read( buffer, BUF_SIZE );
    header = serializer.deserializeStart( BUF_SIZE, buffer );
    resW = header.width;

    if ( header.height % procNum || header.width % procNum )
        throw "Invalid matrix size";
    if ( aLine.width() != header.height )
        throw "Invalid matrices dimensions( non equal )";

    header.height /= procNum;
    myShift = sizeof( SHeader ) + header.height * myID * header.width * ELEMENT_SIZE_BY_TYPE[ header.dataType ];
    bFile.seekg( myShift );
    size = header.height * header.width * ELEMENT_SIZE_BY_TYPE[ header.dataType ];
    matrix<int> bLine( header.height, header.width );
    while ( size > 0 )
    {
        size_t readSize = size < BUF_SIZE ? size : BUF_SIZE;
        bFile.read( buffer, readSize );
        serializer.deserializeStep( buffer, readSize, bLine );
        size -= readSize;
    }
    serializer.deserializeStop();

    aFile.close();
    bFile.close();

    //----------------------------------------------------
 
    double start = MPI_Wtime();

    int blockWidth = bLine.width() / procNum;
    matrix<int> result( aLine.height(), bLine.width() );
    MPI_Datatype SEND_MATRIX_BLOCK;
    MPI_Datatype RECV_MATRIX_BLOCK;
    MPI_Request sendRequest;
    std::vector< MPI_Request > recvRequests;
    std::vector< MPI_Status >  recvStatuses;

    for ( int iter = 0; iter < procNum; ++iter )
    {
        matrix<int> myBlock = getBlock( bLine, procNum, iter );

        if ( iter == 0 )
        {
            checkres( MPI_Type_vector( myBlock.height(), myBlock.width(), myBlock.dataWidth(), MPI_INT, &SEND_MATRIX_BLOCK ) );
            checkres( MPI_Type_commit( &SEND_MATRIX_BLOCK ) );
        }

        for ( int i = 0; i < procNum; ++i )
            if ( i != myID )
                checkres( MPI_Isend( myBlock.shiftedRaw(), 1, SEND_MATRIX_BLOCK, i, BLOCK_BCAST, MPI_COMM_WORLD, &sendRequest) );

        matrix<int> bTranspLine( aLine.width(), blockWidth );
        for ( long i = 0; i < bLine.height(); ++i )
            for ( long q = 0; q < blockWidth; ++q )
                bTranspLine.at( i + myID * bLine.height(), q ) = myBlock.at( i, q );
        
        if ( iter == 0 )
        {
            checkres( MPI_Type_vector( myBlock.height(), myBlock.width(), bTranspLine.dataWidth(), MPI_INT, &RECV_MATRIX_BLOCK ) );
            checkres( MPI_Type_commit( &RECV_MATRIX_BLOCK ) );
        }

        for ( int i = 0; i < procNum; ++i )
        {
            if ( i != myID )
            {
                char* dataPos = (char *)bTranspLine.raw() + i * bLine.height() * blockWidth * ELEMENT_SIZE_BY_TYPE[ bTranspLine.dataType() ];
                recvRequests.push_back( MPI_Request() );
                recvStatuses.push_back( MPI_Status() );
                checkres( MPI_Irecv( dataPos, 1, RECV_MATRIX_BLOCK, i, BLOCK_BCAST, MPI_COMM_WORLD, &recvRequests.back() ) );
            }
        }
        checkres( MPI_Waitall( recvRequests.size(), recvRequests.data(), recvStatuses.data() ) );

        matrix<int>* resBlock = matrix_helper<int>::mul( aLine, bTranspLine );
        long targetCol = iter * blockWidth;
        for ( long i = 0; i < resBlock->height(); ++i )
            for ( long q = 0; q < resBlock->width(); ++q )
                result.at( i, q + targetCol ) = resBlock->at( i, q );    

        delete resBlock;
    }

    double end = MPI_Wtime();
    MASTERPRINT( "TOTAL TIME: " << end - start << "\r\n" );

    //-----------------------------------------------------

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
