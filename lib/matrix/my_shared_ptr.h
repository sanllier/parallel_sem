#ifndef _MY_SHARED_PTR
#define _MY_SHARED_PTR

#include <map>
#include <iterator>

#define nullptr 0

//--------------------------------------------------

template< class T >
class shared_ptr
{
public:
	shared_ptr():m_ptr( nullptr ) {}
	shared_ptr( T* ptr )
	{
		std::map< long long, int >::iterator found = m_map.find( ( long long )ptr );
		m_ptr = ptr;
		if ( found != m_map.end() )
			++( found->second );
		else
			m_map[ ( long long )ptr ] = 1;
	}
	shared_ptr( shared_ptr& ptr )
	{
		std::map< long long, int >::iterator found = m_map.find( ( long long )ptr.m_ptr );
		m_ptr = ptr.m_ptr;
		if ( found != m_map.end() )
			++( found->second );
		else
			m_map[ ( long long )ptr.m_ptr ] = 1;
	}
	~shared_ptr()
	{
		free();
	}

	T* get()
	{
		return m_ptr;
	}

    const T* get() const
    {
        return m_ptr;
    }

	void free()
	{
		std::map< long long, int >::iterator found = m_map.find( ( long long )m_ptr );
		if ( found != m_map.end() )
		{
			--( found->second );
			if ( found->second <= 0 )
            {
                delete (T*)( found->first );
				m_map.erase( found );
            }
		}
		m_ptr = nullptr;
	}

	void reset( T* ptr )
	{
		free();
		
		std::map< long long, int >::iterator found = m_map.find( ( long long )ptr );
		m_ptr = ptr;
		if ( found != m_map.end() )
			++( found->second );
		else
			m_map[ ( long long )ptr ] = 1;
	}

    void operator=( const shared_ptr<T>& ptr )
    {
        free();
        std::map< long long, int >::iterator found = m_map.find( ( long long )ptr.m_ptr );
		m_ptr = ptr.m_ptr;
		if ( found != m_map.end() )
			++( found->second );
		else
			m_map[ ( long long )m_ptr ] = 1;
    }

    T& operator*()
    {
        if ( m_ptr == nullptr )
            throw "Dereferencing NULL pointer";

        return *m_ptr;
    }

private:
	static std::map< long long, int > m_map;
	T* m_ptr;
};

template< class T >
std::map< long long, int > shared_ptr<T>::m_map;

//--------------------------------------------------

#endif
