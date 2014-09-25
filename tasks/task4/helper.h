#ifndef HELPER_H
#define HELPER_H

#include <iostream>
#include <fstream>
#include <vector>

#include "matrix\matrix.h"
#include "matrix\matrix_helper.h"
#include "matrix\matrix_serialization.h"
#include "mpi.h"
#include "parparser\parparser.h"

//----------------------------------------------------------

#define MASTERPRINT( _x_ ) if ( myID == ROOT_ID ) std::cout << _x_ ; std::cout.flush()

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

//----------------------------------------------------------

template< class T >
inline Matrix::matrix<T> getBlock( Matrix::matrix<T>& mat, int procNum, int block )
{
    long blockWidth = mat.width() / procNum;
    matrix<T> subMat( mat.height(), blockWidth );
    subMat.weakSubmatrix( mat, 0, block * blockWidth, mat.height(), blockWidth );
    return subMat;
}

void* deserializeLine( const char* fileName, int lineNum, int totalLines, Matrix::SHeader* fullHeader )
{
    std::ifstream file( fileName, std::ios::binary );
    if ( !file.good() )
    {
        *fullHeader = Matrix::SHeader();
        return 0;
    }

    matrix_serialization serializer;
    size_t size = fileSize( file );
    file.read( buffer, BUF_SIZE );
    Matrix::SHeader header = serializer.deserializeStart( BUF_SIZE, buffer );
    if ( header.height % totalLines || header.width % totalLines )
        throw "Invalid matrix size";

    *fullHeader = header;

    header.height /= totalLines;
    int myShift = sizeof( SHeader ) + header.height * lineNum * header.width * ELEMENT_SIZE_BY_TYPE[ header.dataType ];
    file.seekg( myShift );
    size = header.height * header.width * ELEMENT_SIZE_BY_TYPE[ header.dataType ];

    void* line = matrix_serialization::newMatrixByType( header );
    while ( size > 0 )
    {
        size_t readSize = size < BUF_SIZE ? size : BUF_SIZE;
        file.read( buffer, readSize );
        if ( !matrix_serialization::deserializeStepByType( serializer, header, buffer, readSize, line ) )
            throw "Error while reading the line";

        size -= readSize;        
    }
    serializer.deserializeStop();
    file.close();

    return line;
}


//----------------------------------------------------------

#endif
