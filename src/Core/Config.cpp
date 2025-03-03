#include "Config.h"

#include <fstream>

int _dbglevel = 0;
#ifndef NDEBUG
std::ofstream _tmpos("tmpos.dat", std::ios::binary);
std::ostream tmpos(_tmpos.rdbuf());
#endif  // NDEBUG
