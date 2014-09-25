#ifndef MATRIX_RESOLVER_H
#define MATRIX_RESOLVER_H

#include "matrix.h"
#include "matrix_serialization.h"
#include "matrix_helper.h"

namespace Matrix {
//--------------------------------------------------------------

class matrix_resolver
{
public:
    static matrix_resolver getResolver( const SHeader& header )
    {
        matrix_resolver resolver;
        resolver.m_type = header.dataType;
        return resolver;
    }

    template< typename T > matrix<T>& operator()( void* matr ) { return matrix<int>(); }

public:
    EDataType m_type;
};

template<> matrix<int>& matrix_resolver::operator()( void* matr ) { return *( (matrix<int>*)matr ); }
template<> matrix<long>& matrix_resolver::operator()( void* matr ) { return *( (matrix<long>*)matr ); }
template<> matrix<float>& matrix_resolver::operator()( void* matr ) { return *( (matrix<float>*)matr ); }
template<> matrix<double>& matrix_resolver::operator()( void* matr ) { return *( (matrix<double>*)matr ); }

//--------------------------------------------------------------
}

#endif