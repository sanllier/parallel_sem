#include "parparser.h"

//--------------------------------------------------------------

parparser::parparser( int argc, char** argv )
{
    if ( argv && !parse( argc, argv ) )
        clear();
}

bool parparser::parse( int argc, char** argv )
{
    for ( int i = 1; i < argc; ++i )
    {
        if ( argv[i][0] == '-' )
        {
            if ( argv[i][1] == '-' )
            {
                const char* delimPos = strstr( argv[i], "=" );
                if ( !delimPos )
                    continue;

                std::string key = argv[i] + 2;
                key.resize( delimPos - argv[i] - 2 );
                std::string value = argv[i] - delimPos == 0 ? "" : delimPos + 1;
                
                m_map[ key ] = value;
            }
            else
            {
                if ( i + 1 >= argc )
                    return false;

                std::string key = argv[i] + 1;
                std::string value = argv[ i + 1 ];
                if ( value[0] == '-' )
                    return false;

                ++i;

                m_map[ key ] = value;
            }
        }
    }

    return true;
}

KVPair parparser::get( const char* key )
{
    auto found = m_map.find( key );
    if ( found == m_map.end() )
        return KVPair();
    else
        return KVPair( found->first.c_str(), found->second.c_str() );
}
