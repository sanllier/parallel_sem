#ifndef MATRIX_H
#define MATRIX_H

#include <cstring>
#include <memory>
#include <exception>
#include <complex>

namespace Matrix {
//--------------------------------------------------------------

enum EDataType : int
{
	UNDEFINED_DATA_TYPE,
	INT,
	LONG,
	FLOAT,
	DOUBLE,
    COMPLEX_INT,
    COMPLEX_LONG,
	COMPLEX_FLOAT,
	COMPLEX_DOUBLE
};

enum EMatrixType : int
{
	UNDEFINED_MATRIX_TYPE,
	DENSE_NORMAL
};

//--------------------------------------------------------------

template< class T >
class matrix
{
private:
	struct SMatrixInfo
	{
        long height;
        long width;
		long pivotRow;
		long pivotCol;
		SMatrixInfo( long height = 0, long width = 0, long row = 0, long col = 0)
            : height( height )
            , width( width )
            , pivotRow( row )
            , pivotCol( col ) {}
	};

public:
	matrix()
		: m_data( nullptr )
		, m_dataHeight(0)
        , m_dataWidth(0)
		, m_dataType( getType() )
	    , m_matrixType( DENSE_NORMAL ) {}
	matrix( long height, long width )		
        : m_info( height, width ) 
        , m_dataHeight( height )
        , m_dataWidth( width )
		, m_dataType( getType() )
		, m_matrixType( DENSE_NORMAL )
	{
		if ( m_dataHeight && m_dataWidth )
		{
			m_data.reset( new T[ m_dataHeight * m_dataWidth ] );	
            clear();
		}
		else
		{
			m_info = SMatrixInfo();
            m_dataHeight = 0;
            m_dataWidth  = 0;           
		}
	}
	matrix( const matrix& mat )
        : m_info( mat.height(), mat.width() )
		, m_dataHeight( mat.m_info.height )
        , m_dataWidth( mat.m_info.width )
		, m_dataType( getType() )
		, m_matrixType( DENSE_NORMAL )
	{
        if ( m_dataHeight && m_dataWidth )
		{
			m_data.reset( new T[ m_dataHeight * m_dataWidth ] );	
            for ( long i = 0; i < mat.height(); ++i )
            {
                const size_t thisShift = m_dataWidth * i;
                const size_t matShift  = mat.m_dataWidth * ( i + mat.m_info.pivotRow ) + mat.m_info.pivotCol;
                std::memcpy( m_data.get() + thisShift, mat.m_data.get() + matShift, m_dataWidth * sizeof(T) );
            }            
		}	
		else
		{
			m_info = SMatrixInfo();
            m_dataHeight = 0;
            m_dataWidth  = 0;           
		}
	}
    ~matrix() {}

	//----------------------- COPY ------------------------
	void weakCopy( const matrix& mat )
	{
        m_info       = mat.m_info;
        m_dataHeight = mat.m_dataHeight;
        m_dataWidth  = mat.m_dataWidth;
        m_dataType   = mat.m_dataType;
        m_matrixType = mat.m_matrixType;
		m_data       = mat.m_data;
	}
	void strongCopy( const matrix& mat )
	{
        matrix<T> temp( mat );
        this->weakCopy( temp );
	}
    void equalStrongCopy( const matrix& mat )
    {
        if ( this->height() != mat.height() || this->width() != mat.width() )
            throw "equalStrongCopy error";

        // CRAP
        for ( long i = 0; i < mat.height(); ++i )
            for ( long q = 0; q < mat.width(); ++q )
            this->at( i, q ) = mat.at( i, q );
    }
	matrix& operator=( const matrix& mat )
	{
		weakCopy( mat );
		return *this;
	}

	//---------------------- ACCESS -----------------------
	inline T& at( long row, long col )
	{
        #ifndef _PERFORMANCE_ 
		    if ( row >= height() )
			    throw std::out_of_range( "matrix row out of range" );
		    else if ( col >= width() )
			    throw std::out_of_range( "matrix column out of range" );
        #endif

        const size_t shift = m_dataWidth * ( row + m_info.pivotRow ) + ( col + m_info.pivotCol );
		return m_data.get()[ shift ];
	}
    inline const T& at( long row, long col ) const
	{
        #ifndef _PERFORMANCE_ 
		if ( row >= height() )
			throw std::out_of_range( "matrix row out of range" );
		else if ( col >= width() )
			throw std::out_of_range( "matrix column out of range" );
        #endif

        const size_t shift = m_dataWidth * ( row + m_info.pivotRow ) + ( col + m_info.pivotCol );
		return m_data.get()[ shift ];
	}
    inline long height() const { return m_info.height; }
    inline long width()  const { return m_info.width; }
    inline long dataHeight() const { return m_dataHeight; }
    inline long dataWidth()  const { return m_dataWidth; }
    inline EDataType dataType() const { return m_dataType; }
    inline void setMatrixType( EMatrixType type ) { m_matrixType = type; }
    inline EMatrixType matrixType() const { return m_matrixType; }
    inline void* raw() { return m_data.get(); }
    inline const void* raw() const { return m_data.get(); }
    inline void* shiftedRaw() { return m_data.get() + m_info.pivotRow * dataWidth() + m_info.pivotCol; } 
    inline const void* shiftedRaw() const { return m_data.get() + m_info.pivotRow * dataWidth() + m_info.pivotCol; } 

    //-----------------------------------------------------
    matrix<T>& strongSubmatrix( const matrix<T>& mat, long row, long col, long height, long width )
    {
        if ( row + height > mat.height() || col + width > mat.width() )
            return *this;

        matrix<T> weakMatr;        
        matrix<T> temp( weakMatr.weakSubmatrix( mat, row, col, height, width ) );
        this->strongCopy( temp );
        return *this;
    }
    matrix<T>& weakSubmatrix( const matrix<T>& mat, long row, long col, long height, long width )
    {
        if ( row + height > mat.height() || col + width > mat.width() )
            return *this;

        m_info = SMatrixInfo( height, width, row, col );
        m_dataHeight = mat.m_dataHeight;
        m_dataWidth  = mat.m_dataWidth;
        m_matrixType = mat.m_matrixType;
        m_data = mat.m_data;
        return *this;
    }
    matrix<T>& insertSubmatrix( const matrix<T>& mat, long row, long col )
    {
        const long insH = row + mat.height() > this->height() ? mat.height() - row : row + mat.height();
        const long insW = col + mat.width()  > this->width() ? mat.width() - col : col + mat.width();

        for ( long i = row; i < insH; ++i )
            for ( long q = col; q < insW; ++q )
                this->at( i, q ) = mat.at( i - row, q - col );

        return *this;
    }
    void clear()
    {
        std::memset( m_data.get(), 0, m_dataHeight * m_dataWidth * sizeof(T) );
    }

    //-------------------- OPERATIONS ---------------------
    matrix<T>& add( const matrix<T>& matr )
    {
        const long height = this->height();
        const long width  = this->width();

        if ( height != matr.height() || width != matr.width() )
            return *this;

        for ( long i = 0; i < height; ++i )
            for ( long q = 0; q < width; ++q )
                at( i, q ) += matr.at( i, q );

        return *this;
    }
    matrix<T>& mul( const matrix<T>& matr )
    {
        const long height = this->height();
        const long width  = this->width();
        if ( width != matr.height() )
            return *this;

        matrix<T> temp;
        temp.strongCopy( *this );
        m_data.reset( new T[ height * matr.width() ] );
        m_info = SMatrixInfo( height, matr.width() );
        m_dataHeight = height;
        m_dataWidth  = matr.width();

        T acc = T();
        for ( long i = 0; i < height; ++i )
        {
            for ( long q = 0; q < matr.width(); ++q )
            {
                for ( long k = 0; k < width; ++k )
                    acc += temp.at( i, k ) * matr.at( k, q );
 
                at( i, q ) = acc;
                acc = T();
            }
        }
        return *this;
    }

private:
	EDataType getType() { return UNDEFINED_DATA_TYPE; }

private:
    SMatrixInfo m_info;
	long m_dataHeight;
	long m_dataWidth;
	EDataType   m_dataType;
	EMatrixType m_matrixType;

	std::shared_ptr<T> m_data;
};

template<> EDataType matrix<int>::getType() { return INT; }
template<> EDataType matrix<long>::getType() { return LONG; }
template<> EDataType matrix<float>::getType() { return FLOAT; }
template<> EDataType matrix<double>::getType() { return DOUBLE; }
template<> EDataType matrix< std::complex<int> >::getType() { return COMPLEX_INT; }
template<> EDataType matrix< std::complex<long> >::getType() { return COMPLEX_LONG; }
template<> EDataType matrix< std::complex<float> >::getType() { return COMPLEX_FLOAT; }
template<> EDataType matrix< std::complex<double> >::getType() { return COMPLEX_DOUBLE; }

//--------------------------------------------------------------
}

#endif
