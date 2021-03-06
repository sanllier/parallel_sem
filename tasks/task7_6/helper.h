#ifndef HELPER_H
#define HELPER_H

#include <iostream>
#include <fstream>
#include <vector>

#include "matrix/matrix.h"
#include "matrix/matrix_helper.h"
#include "matrix/matrix_serialization.h"
#include "mpi.h"
#include "parparser/parparser.h"

//----------------------------------------------------------

#define MASTERPRINT( _x_ ) if ( myID == ROOT_ID ) std::cout << _x_ ; std::cout.flush()

#define RESOLVE_OP( _header_, _mat_, _op_ )     \
    if ( _header_.dataType == INT )             \
        ( (matrix<int>*)_mat_ )->_op_;          \
    else if ( _header_.dataType == LONG )       \
        ( (matrix<long>*)_mat_ )->_op_;         \
    else if ( _header_.dataType == FLOAT )      \
        ( (matrix<float>*)_mat_ )->_op_;        \
    else if ( _header_.dataType == DOUBLE )     \
        ( (matrix<double>*)_mat_ )->_op_;       \

#define RESOLVE_VAR( _header_, _mat_, _op_, _var_ ) \
    if ( _header_.dataType == INT )                 \
        _var_ = ( (matrix<int>*)_mat_ )->_op_;      \
    else if ( _header_.dataType == LONG )           \
        _var_ = ( (matrix<long>*)_mat_ )->_op_;     \
    else if ( _header_.dataType == FLOAT )          \
        _var_ = ( (matrix<float>*)_mat_ )->_op_;    \
    else if ( _header_.dataType == DOUBLE )         \
        _var_ = ( (matrix<double>*)_mat_ )->_op_;   \

//----------------------------------------------------------

enum
{
    ROOT_ID = 0,
    BLOCK_BCAST
};

static const size_t BUF_SIZE = sizeof( Matrix::SHeader );
static char buffer[ BUF_SIZE ];

//----------------------------------------------------------

inline void checkres( int res )
{
    if ( res != MPI_SUCCESS )
    {
        char estring[ MPI_MAX_ERROR_STRING ];
        int len;
        MPI_Error_string( res, estring, &len );
        throw estring;
    }
}

inline MPI_Datatype MPI_Type_vector_wrapper( int count, int blocklength, int stride, MPI_Datatype oldtype )
{
    MPI_Datatype temp;
    checkres( MPI_Type_vector( count, blocklength, stride, oldtype, &temp ) );
    checkres( MPI_Type_commit( &temp ) );
    return temp;
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
bool parallelMatrixSerialization( const Matrix::SHeader& fullHeader, const Matrix::SHeader& partHeader, const Matrix::matrix<T>& mat, const char* fileName, MPI_Comm comm )
{
    if ( !fileName || !fileName[0] )
        return false;

    int procNum = 0;
    int myID = 0;
    checkres( MPI_Comm_size( comm, &procNum ) );
    checkres( MPI_Comm_rank( comm, &myID ) );

    for ( int i = 0; i < procNum; ++i )
    {
        if ( myID == i )
        {
            if ( myID == ROOT_ID )
            {
                std::ofstream res( fileName, std::ios::binary );
                if ( !res.good() )
                    return false;

                res.write( (const char*)&fullHeader, sizeof( Matrix::SHeader ) );
                const char* raw = (char*)mat.raw();
                res.write( raw, ELEMENT_SIZE_BY_TYPE[ partHeader.dataType ] * mat.height() * mat.width() );
                res.close();
            }
            else
            {
                std::fstream res;
                res.open ( fileName, std::fstream::out | std::fstream::app | std::ios::binary );
                if ( !res.good() )
                    return false;

                const char* raw = (char*)mat.raw();
                res.write( raw, ELEMENT_SIZE_BY_TYPE[ partHeader.dataType ] * mat.height() * mat.width() );
                res.close();
            }
        }
        MPI_Barrier( comm );
    }

    return true;
};

//----------------------------------------------------------

template< class T >
inline Matrix::matrix<T> getBlock( Matrix::matrix<T>& mat, int procNum, int block )
{
    long blockWidth = mat.width() / procNum;
    Matrix::matrix<T> subMat( mat.height(), blockWidth );
    subMat.weakSubmatrix( mat, 0, block * blockWidth, mat.height(), blockWidth );
    return subMat;
}

template< class T >
Matrix::matrix<T> deserializeLine( const char* fileName, int lineNum, int totalLines, Matrix::SHeader* fullHeader, Matrix::SHeader* lineHeader )
{
    std::ifstream file( fileName, std::ios::binary );
    if ( !file.good() )
    {
        *fullHeader = Matrix::SHeader();
        *lineHeader = Matrix::SHeader();
        return Matrix::matrix<T>();
    }

    Matrix::matrix_serialization serializer;
    size_t size = fileSize( file );
    file.read( buffer, BUF_SIZE );
    Matrix::SHeader header = serializer.deserializeStart( BUF_SIZE, buffer );
    if ( header.height % totalLines || header.width % totalLines )
        throw "Invalid matrix size";

    *fullHeader = header;

    header.height /= totalLines;
    int myShift = sizeof( Matrix::SHeader ) + header.height * lineNum * header.width * ELEMENT_SIZE_BY_TYPE[ header.dataType ];
    file.seekg( myShift );
    size = header.height * header.width * ELEMENT_SIZE_BY_TYPE[ header.dataType ];

    *lineHeader = header;

    Matrix::matrix<T> line( header.height, header.width );
    while ( size > 0 )
    {
        size_t readSize = size < BUF_SIZE ? size : BUF_SIZE;
        file.read( buffer, readSize );
        serializer.deserializeStep<T>( buffer, readSize, line );
        size -= readSize;        
    }
    serializer.deserializeStop();
    file.close();

    return line;
}

//----------------------------------------------------------

template< class T, MPI_Datatype MPI_TYPE >
Matrix::matrix<T> parallelRead( const char* file, int rank, int xDim, int yDim, const int coords[2], MPI_Comm comm, Matrix::SHeader* header )
{
    Matrix::matrix<T> block;
    if ( !file || !file[0] )
        return block;

    Matrix::SHeader tempHeader;
    if ( rank == ROOT_ID )
    {
        std::ifstream inpFile( file, std::ios::binary );
        if ( !inpFile.good() )
            return block;

        inpFile.read( (char *)( &tempHeader ), sizeof( tempHeader ) );
        inpFile.close();

        checkres( MPI_Bcast( &tempHeader, sizeof( tempHeader ), MPI_CHAR, ROOT_ID, MPI_COMM_WORLD ) );
    }
    else
    {
        checkres( MPI_Bcast( &tempHeader, sizeof( tempHeader ), MPI_CHAR, ROOT_ID, MPI_COMM_WORLD ) );
    }

    if ( header )
        *header = tempHeader;

    block = Matrix::matrix<T>( tempHeader.height / xDim, tempHeader.width / yDim );
    MPI_Datatype blockType = MPI_Type_vector_wrapper( block.height(), block.width(), tempHeader.width, MPI_TYPE );
    MPI_File matFile;
    MPI_Status status;
    checkres( MPI_File_open( comm, const_cast< char* >(file), MPI_MODE_RDONLY,  MPI_INFO_NULL, &matFile ) );

    MPI_Offset myOffset = ( coords[0] * block.height() * tempHeader.width + coords[1] * block.width() ) * sizeof( T );
    checkres( MPI_File_set_view( matFile, sizeof( Matrix::SHeader ) + myOffset , MPI_TYPE, blockType, "native", MPI_INFO_NULL) );
    MPI_File_read_all( matFile, block.raw(), block.height() * block.width(), MPI_TYPE, &status );
    MPI_File_close( &matFile );

    return block;
}

template< class T, MPI_Datatype MPI_TYPE >
bool parallelWrite( const  char* file, int rank, int xDim, int yDim, const int coords[2], MPI_Comm comm, const Matrix::matrix<T>& block )
{
    if ( !file || !file[0] )
        return false;

    Matrix::SHeader header = { block.height() * xDim, block.width() * yDim, block.dataType(), block.matrixType() };
    if ( rank == ROOT_ID )
    {
        std::ofstream outFile( file, std::ios::binary );
        if ( !outFile.good() )
            return false;

        outFile.write( ( char* )( &header ), sizeof( header ) );
        outFile.close();
    }
    checkres( MPI_Barrier( comm ) );

    MPI_Datatype blockType = MPI_Type_vector_wrapper( block.height(), block.width(), header.width, MPI_TYPE );
    MPI_File matFile;
    MPI_Status status;
    checkres( MPI_File_open( comm, const_cast< char* >(file), MPI_MODE_CREATE | MPI_MODE_WRONLY,  MPI_INFO_NULL, &matFile ) );

    MPI_Offset myOffset = ( coords[0] * block.height() * header.width + coords[1] * block.width() ) * sizeof( T );
    checkres( MPI_File_set_view( matFile, sizeof( Matrix::SHeader ) + myOffset, MPI_TYPE, blockType, "native", MPI_INFO_NULL) );
    MPI_File_write_all( matFile, const_cast< void* >( block.raw() ), block.height() * block.width(), MPI_TYPE, &status );
    MPI_File_close( &matFile );

    return true;
}

#endif
