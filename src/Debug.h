#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>
#ifdef DEBUG
#define DEBUG_MSG(str) do { std::cout << "DEBUG: " << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

#endif
