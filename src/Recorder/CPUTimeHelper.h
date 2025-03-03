
#include <chrono>

#pragma once

struct CPUTimeHelper {
  using HRC = std::chrono::high_resolution_clock;
  std::chrono::time_point<HRC> start;
  CPUTimeHelper() { reset(); }
  void reset() { start = HRC::now(); }
  double operator()() const {
    std::chrono::duration<double> e = HRC::now() - start;
    return e.count();
  }
};