#ifndef CONFIG_H
#define CONFIG_H

#include <unistd.h>

#include <iostream>

#include "Recorder/SpdlogWrapper.h"

using Real = double;

#ifndef ORDER
#define ORDER 4
#endif

#ifndef DIM
#define DIM 3
#endif

#if DIM == 2
const int SpaceDim = 2;
#elif DIM == 3
const int SpaceDim = 3;
#endif
extern int _dbglevel;

//=================================================

// debug output issues
#ifndef NDEBUG

extern std::ostream tmpos;
#define push_dbglevel(x) _dbglevel = (_dbglevel << 4) | (x)
#define pop_dbglevel() _dbglevel = _dbglevel >> 4
#define reset_dbglevel(x) _dbglevel = (x)
#define get_dbglevel() (_dbglevel & 0x0f)

#else

#define tmpos 0 && std::cout
#define push_dbglevel(x)
#define pop_dbglevel()
#define reset_dbglevel(x)
#ifndef DBGLEVEL
#define DBGLEVEL 1
#endif
#define get_dbglevel() (DBGLEVEL)

#endif  // NDEBUG

#define dbgcout (get_dbglevel() >= 0) && std::cout
#define dbgcout1 (get_dbglevel() >= 1) && std::cout << "  "
#define dbgcout2 (get_dbglevel() >= 2) && std::cout << "    "
#define dbgcout3 (get_dbglevel() >= 3) && std::cout << "      "

inline void breakpoint() { __asm__ __volatile__("int3"); }

inline void keepSleep() {
  InfoMessage("Start Sleep ...");
  volatile bool stop = true;
  while (stop) {
    sleep(1);
  }
}

inline Real distTol(Real tol = -1) noexcept {
  static Real distTol = 1e-12;
  if (tol > 0) distTol = tol;
  return distTol;
}

inline Real newtonTol(Real tol = -1) noexcept {
  static Real newtonTol = 1e-15;
  if (tol > 0) newtonTol = tol;
  return newtonTol;
}


inline Real newtonSubCount(int Count = -1) noexcept {
  static int newtonSubCount = 200;
  if (Count > 0) newtonSubCount = Count;
  return newtonSubCount;
}

inline int newtonMaxIter(int num = -1) noexcept {
  static int newtonMaxIter = 20;
  if (num > 0) newtonMaxIter = num;
  return newtonMaxIter;
}

inline int defaultHighPrecision = 128;

#endif  // CONFIG_H
