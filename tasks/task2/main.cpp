#include <iostream>
#include <fstream>
#include <ctime>

#include "matrix\matrix.h"
#include "matrix\matrix_helper.h"
#include "matrix\matrix_serialization.h"
//#include "matrix\matrix_mpi.h"
#include "mpi.h"
#include "parparser\parparser.h"

//----------------------------------------------------------

using namespace Matrix;

//----------------------------------------------------------

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
    int* blockBuf = new int[ bLine.height() * blockWidth ];
    matrix<int> result( aLine.height(), bLine.width() );

    for ( int iter = 0; iter < procNum; ++iter )
    {
        matrix<int> myBlock = getBlock( bLine, procNum, iter );
        long pos = 0;
        for ( long i = 0; i < myBlock.height(); ++i )
            for ( long q = 0; q < myBlock.width(); ++q )
                blockBuf[ pos++ ] = myBlock.at( i, q );

        for ( int i = 0; i < procNum; ++i )
        {
            if ( i == myID )
                continue;
            checkres( MPI_Send( blockBuf, pos, MPI_INT, i, BLOCK_BCAST, MPI_COMM_WORLD ) );
        }

        matrix<int> bTranspLine( aLine.width(), blockWidth );
        for ( long i = 0; i < bLine.height(); ++i )
            for ( long q = 0; q < blockWidth; ++q )
                bTranspLine.at( i + myID * bLine.height(), q ) = myBlock.at( i, q );

        // FIX
        //int recvCount = procNum - 1;
        //while ( recvCount > 0 )
        //{
        //    MPI_Status recvStatus;
        //    MPI_Recv( blockBuf, pos, MPI_INT, MPI_ANY_SOURCE, BLOCK_BCAST, MPI_COMM_WORLD, &recvStatus);

        //    long temp = 0;
        //    long targetRow = recvStatus.MPI_SOURCE * bLine.height();
        //    for ( long i = 0; i < bLine.height(); ++i )
        //        for ( long q = 0; q < blockWidth; ++q )
        //            bTranspLine.at( i + targetRow, q ) = blockBuf[ temp++ ];

        //    --recvCount;
        //}
        // CRAP
        for ( int q = 0; q < procNum; ++q )
        {
            if ( q == myID )
                continue;

            MPI_Status recvStatus;
            MPI_Recv( blockBuf, pos, MPI_INT, q, BLOCK_BCAST, MPI_COMM_WORLD, &recvStatus);

            long temp = 0;
            long targetRow = q * bLine.height();
            for ( long i = 0; i < bLine.height(); ++i )
                for ( long k = 0; k < blockWidth; ++k )
                    bTranspLine.at( i + targetRow, k ) = blockBuf[ temp++ ];

        }

        matrix<int>* resBlock = matrix_helper<int>::mul( aLine, bTranspLine );
        long targetCol = iter * blockWidth;
        for ( long i = 0; i < resBlock->height(); ++i )
            for ( long q = 0; q < resBlock->width(); ++q )
                result.at( i, q + targetCol ) = resBlock->at( i, q );    

        delete resBlock;
    }

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

    double end = MPI_Wtime();
    if ( myID == ROOT_ID )
        std::cout << "TOTAL TIME: " << end - start << "\r\n";

    //-----------------------------------------------------
    
    checkres( MPI_Finalize() );
    return 0;
}
