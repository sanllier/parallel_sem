#ifndef MATRIX_SERIALIZATION_H
#define MATRIX_SERIALIZATION_H

#include "matrix.h"

#include <fstream>

//-------------------------------------------------------------

#ifdef _MSC_VER 
    #ifdef _WIN64
        static const int ELEMENT_SIZE_BY_TYPE[] = { 0, 4, 8, 4, 8, 8, 16, 8, 16 };
    #else
        static const int ELEMENT_SIZE_BY_TYPE[] = { 0, 4, 4, 4, 8, 8, 8, 8, 16 };
    #endif
#endif

#ifdef __GNUC__
    #ifdef __x86_64__
        static const int ELEMENT_SIZE_BY_TYPE[] = { 0, 4, 8, 4, 8, 8, 16, 8, 16 };
    #else
        static const int ELEMENT_SIZE_BY_TYPE[] = { 0, 4, 4, 4, 8, 8, 8, 8, 16 };
    #endif
#endif

//-------------------------------------------------------------

namespace Matrix {
//--------------------------------------------------------------

struct SHeader
{
    long height;
	long width;
	EDataType   dataType;
	EMatrixType matrixType;
};

//--------------------------------------------------------------

class matrix_serialization
{
public:
    matrix_serialization():m_buf(0), m_bufSize(0), m_dataSize(0), m_pos(0) {}
    ~matrix_serialization()
    {
        if ( m_buf )
        {
            delete[] m_buf;
            m_buf = 0;
        }
    }

    //---------------------- HELPER -----------------------

    static void* newMatrixByType( const SHeader& header )
    {
        void* temp = 0;
        switch ( header.dataType )
		{
		case INT:
			temp = new matrix<int>( header.height, header.width );
			break;
		case LONG:
			temp = new matrix<long>( header.height, header.width );
			break;
		case FLOAT:
			temp = new matrix<float>( header.height, header.width );
			break;
		case DOUBLE:
			temp = new matrix<double>( header.height, header.width );
			break;
		case COMPLEX_INT:
			temp = new matrix< std::complex<int> >( header.height, header.width );
			break;
		case COMPLEX_LONG:
			temp = new matrix< std::complex<long> >( header.height, header.width );
			break;
		case COMPLEX_FLOAT:
			temp = new matrix< std::complex<float> >( header.height, header.width );
			break;
		case COMPLEX_DOUBLE:
			temp = new matrix< std::complex<double> >( header.height, header.width );
			break;
		case UNDEFINED_DATA_TYPE:
			return 0;
		default:
			return 0;
		}
    }

    static bool deserializeStepByType( matrix_serialization& serializer, const SHeader& header, const void* buf, size_t& size, void* mat )
    {
        if ( !mat )
            return false;

		switch ( header.dataType )
		{
		case INT:
            serializer.deserializeStep( buf, size, *( (matrix<int>*)mat ) ); 
			break;
		case LONG:
			serializer.deserializeStep( buf, size, *( (matrix<long>*)mat ) ); 
			break;
		case FLOAT:
			serializer.deserializeStep( buf, size, *( (matrix<float>*)mat ) ); 
			break;
		case DOUBLE:
			serializer.deserializeStep( buf, size, *( (matrix<double>*)mat ) ); 
			break;
		case COMPLEX_INT:
			serializer.deserializeStep( buf, size, *( (matrix< std::complex<int> >*)mat ) ); 
			break;
		case COMPLEX_LONG:
			serializer.deserializeStep( buf, size, *( (matrix< std::complex<long> >*)mat ) ); 
			break;
		case COMPLEX_FLOAT:
			serializer.deserializeStep( buf, size, *( (matrix< std::complex<float> >*)mat ) ); 
			break;
		case COMPLEX_DOUBLE:
			serializer.deserializeStep( buf, size, *( (matrix< std::complex<double> >*)mat ) ); 
			break;
		case UNDEFINED_DATA_TYPE:
			return false;			
		default:
			return false;	
		}

        return true;
    }

    //------------------- SERIALIZATION -------------------
    template< class T >
	const void* seriaizeStart( matrix<T>& mat, size_t bufSize )
	{
		if ( bufSize > 0 )
		{
			if ( m_buf )
				delete[] m_buf;

            SHeader header = { mat.height(), mat.width(), mat.dataType(), mat.matrixType() };
			m_bufSize = bufSize;
			m_buf = new char[ m_bufSize ];
            std::memcpy( m_buf, (char*)&header, sizeof( SHeader ) ); 
			m_dataSize = mat.height() * mat.width() * sizeof(T);
			m_pos = 0;
			return m_buf;
		}
		return 0;
	}
    template< class T >
	const void* serializeStep( size_t& size, matrix<T>& mat )
	{
		if ( !m_buf || !m_bufSize )
		{
			size = 0;
			return 0;
		}

		if ( m_pos < m_dataSize )
		{
			size_t temp = m_bufSize;
			if ( m_pos + temp > m_dataSize )
				temp -= ( m_pos + temp ) - m_dataSize;

			std::memcpy( m_buf, (char*)mat.raw() + m_pos, temp );
			m_pos += temp;
            size = temp;
		}
        else
        {
            size = 0;
            std::memset( m_buf, 0, m_bufSize );
        }
		return m_buf;				
	}
	void serializeStop()
	{
		if ( m_buf )
		{
			delete[] m_buf;
			m_buf = 0;
			m_bufSize  = 0;
            m_dataSize = 0;
			m_pos = 0;
		}
	}

    //------------------ DESERIALIZATION ------------------
	SHeader deserializeStart( size_t bufSize, const void* data )
	{
        if ( bufSize < sizeof( SHeader ) )
            return SHeader();

        SHeader header;
        std::memcpy( &header, data, sizeof( SHeader ) );
        m_dataSize = header.height * header.width * ELEMENT_SIZE_BY_TYPE[ header.dataType ];

		if ( m_buf )
		    delete[] m_buf;

        m_bufSize = bufSize;
		m_buf = new char[ bufSize ];
		m_pos = 0;
        return header;
	}
	template< class T >
    void deserializeStep( const void* data, size_t& dataSize, matrix<T>& mat )
	{
		if ( data )
		{
            size_t temp = m_bufSize > m_dataSize - m_pos ? m_dataSize - m_pos : m_bufSize;
            if ( temp > dataSize ) temp = dataSize;
			std::memcpy( (char*)mat.raw() + m_pos, data, temp );
			m_pos += temp;
			dataSize = temp;
		}
	}
	void deserializeStop()
	{
		if ( m_buf )
		{
			delete[] m_buf;
			m_buf = 0;
			m_bufSize  = 0;
            m_dataSize = 0;
			m_pos = 0;
		}
	}

    //--------------------- HELPERS -----------------------
	template< class T >
	bool writeBinary( const char* fileName, matrix<T>& matr )
	{
		static const size_t BUF_SIZE = sizeof( SHeader ); // CRAP?

		if ( !fileName || !fileName[0] )
			return false;

		std::ofstream oFstr( fileName, std::ios::binary );
		if ( !oFstr.good() )
			return false;

        const void* buf = 0;
		buf = seriaizeStart( matr, BUF_SIZE );
        oFstr.write( (const char*)buf, sizeof( SHeader ) );
		size_t size = 0;		
		do
		{			
			buf = serializeStep( size, matr );
            if ( size > 0 )
			    oFstr.write( (const char*)buf, size );
		} while( size );
		serializeStop();

		oFstr.close();
		return true;
	}

	void* readBinary( const char* fileName, SHeader& header )
	{
		static const size_t BUF_SIZE = sizeof( SHeader ); // CRAP?

		if ( !fileName || !fileName[0] )
			return 0;

		std::ifstream iFstr( fileName, std::ios::binary );
		if ( !iFstr.good() )
			return 0;

		iFstr.seekg (0, iFstr.end);
		size_t length = (size_t)iFstr.tellg();
		iFstr.seekg (0, iFstr.beg);
		if ( length < sizeof( SHeader ) )
			return 0;

		char* buf = new char[ BUF_SIZE ];
		if ( !buf )
			return 0;

		void* undefMatr;

		iFstr.read( buf, BUF_SIZE );
		length -= BUF_SIZE;
		header = deserializeStart( BUF_SIZE, buf );

        undefMatr = newMatrixByType( header );
        if ( !undefMatr )
        { 
			delete[] buf; 
			iFstr.close(); 
			deserializeStop();
        }

		size_t size = 0;
		do
		{
			size = BUF_SIZE > length ? length : BUF_SIZE;
			iFstr.read( buf, size );
			length -= BUF_SIZE;

            if ( !matrix_serialization::deserializeStepByType( *this, header, buf, size, undefMatr ) )
            {
                delete undefMatr;
				undefMatr = 0;
				delete[] buf; 
				iFstr.close(); 
				deserializeStop();
				return 0;	
            }			
		} while( size );
		deserializeStop();

		delete[] buf; 
		iFstr.close(); 
	    deserializeStop();
		return undefMatr;
	}

private:
    char* m_buf;
    size_t m_bufSize;
    size_t m_dataSize;
    size_t m_pos;
};

//--------------------------------------------------------------
}

#endif
