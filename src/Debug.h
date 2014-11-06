/* Any copyright is dedicated to the Public Domain.
 *  * http://creativecommons.org/publicdomain/zero/1.0/ */

#ifndef DEBUG_H
#define DEBUG_H

#include <iostream>
#include <iomanip>

#ifdef DEBUG
#define DEBUG_MSG(str) do { std::cout << std::setprecision(12) << "DEBUG: " << str << std::endl; } while( false )
#else
#define DEBUG_MSG(str) do { } while ( false )
#endif

#endif
