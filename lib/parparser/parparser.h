#ifndef MATRIX_MPI_H
#define MATRIX_MPI_H

#include <map>
#include <string>

//--------------------------------------------------------------

struct KVPair
{
    const char* key;
    const char* value;
    bool isSet;

    KVPair( const char* key = 0, const char* value = 0 ):key(0), value(0), isSet( false ) 
    {
        if ( key && value )
            set( key, value );
    }
    inline void set( const char* key, const char* value ) 
    { 
        this->key = key; 
        this->value = value; 
        isSet = true; 
    }

    bool asBool( bool def = false ) { return isSet ? ( 0 == strcmp( "t", value ) || 0 == strcmp( "true", value ) ) : def; }
    int asInt( int def = 0 ) { return isSet ? atoi( value ) : def; }
    long asLong( long def = 0 ) { return isSet ? atol( value ) : def; }
    float asFloat( float def = 0.0f ) { return isSet ? (float)atof( value ) : def; }
    double asDouble( double def = 0.0 ) { return isSet ? atof( value ) : def; }

    bool check( const char* pattern ) { return ( isSet && pattern ) ? 0 == strcmp( pattern, value ) : false; }
};

class parparser
{
public:
    parparser() {}
    parparser( int argc, char** argv );
    ~parparser() {}

    bool parse( int argc, char** argv );
    void clear() { m_map.clear(); }

    KVPair get( const char* key );

private:
    std::map< std::string, std::string > m_map;
};

#endif
