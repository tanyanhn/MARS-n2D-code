#ifndef _RKBUTCHER_H_
#define _RKBUTCHER_H_

#include <cmath>

#include "Core/Config.h"

namespace RK {
enum Type_Major { ERK = 0, DIRK, ARK, nRK_Family };

enum Type_Minor {
  ForwardEuler = 0,
  ClassicRK4,
  Verner6,
  PrinceDormand8,
  SDIRK2,
  ESDIRK4,
  nRK_Type
};
}  // namespace RK

template <RK::Type_Major Type1, RK::Type_Minor Type2>
struct ButcherTableau;

template <>
struct ButcherTableau<RK::ERK, RK::ForwardEuler> {
  static constexpr int order = 1;
  static constexpr int nStages = 1;
  static constexpr long double a[][1] = {{0}};
  static constexpr long double b[] = {1};
  static constexpr long double c[] = {0};
};

template <>
struct ButcherTableau<RK::ERK, RK::ClassicRK4> {
  static constexpr int order = 4;
  static constexpr int nStages = 4;
  static constexpr long double a[][4] = {
      {0, 0, 0, 0}, {0.5, 0, 0, 0}, {0, 0.5, 0, 0}, {0, 0, 1, 0}};
  static constexpr long double b[] = {1.0 / 6, 1.0 / 3, 1.0 / 3, 1.0 / 6};
  static constexpr long double c[] = {0.0, 0.5, 0.5, 1.0};
};

template <>
struct ButcherTableau<RK::ERK, RK::Verner6> {
  static constexpr int order = 6;
  static constexpr int nStages = 8;

  static constexpr long double a[][8] = {
      {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},                     // Stage 1
      {1.0 / 6.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},               // Stage 2
      {4.0 / 75.0, 16.0 / 75.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},      // Stage 3
      {5.0 / 6.0, -8.0 / 3.0, 5.0 / 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},  // Stage 4
      {-165.0 / 64.0, 55.0 / 6.0, -425.0 / 64.0, 85.0 / 96.0, 0.0, 0.0, 0.0,
       0.0},  // Stage 5
      {12.0 / 5.0, -8.0, 4015.0 / 612.0, -11.0 / 36.0, 88.0 / 255.0, 0.0, 0.0,
       0.0},  // Stage 6
      {-8263.0 / 15000.0, 124.0 / 75.0, -643.0 / 680.0, -81.0 / 250.0,
       2484.0 / 10625.0, 0.0, 0.0, 0.0},  // Stage 7
      {3501.0 / 1720.0, -300.0 / 43.0, 297275.0 / 52632.0, -319.0 / 2322.0,
       24068.0 / 84065.0, 0.0, 3850.0 / 26703.0, 0.0}  // Stage 8
  };

  static constexpr long double b[] = {
      3.0 / 40.0,     0.0, 875.0 / 2244.0,  23.0 / 72.0,
      264.0 / 1955.0, 0.0, 125.0 / 11592.0, 43.0 / 616.0};

  static constexpr long double c[] = {0.0,       1.0 / 6.0, 4.0 / 15.0, 2.0 / 3.0,
                               5.0 / 6.0, 1.0,       1.0 / 15.0, 1.0};
};

template <>
struct ButcherTableau<RK::ERK, RK::PrinceDormand8> {
  static constexpr int order = 8;
  static constexpr int nStages = 13;

  static constexpr long double a[][13] = {
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1.0 / 18, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1.0 / 48, 1.0 / 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {1.0 / 32, 0, 3.0 / 32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {5.0 / 16, 0, -75.0 / 64, 75.0 / 64, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      {3.0 / 80, 0, 0, 3.0 / 16, 3.0 / 20, 0, 0, 0, 0, 0, 0, 0, 0},
      {29443841.0 / 614563906, 0, 0, 77736538.0 / 692538347,
       -28693883.0 / 1125000000, 23124283.0 / 1800000000, 0, 0, 0, 0, 0, 0, 0},
      {16016141.0 / 946692911, 0, 0, 61564180.0 / 158732637,
       22789713.0 / 633445777, 545815736.0 / 2771057229,
       -180193667.0 / 1043307555, 0, 0, 0, 0, 0, 0},
      {39632708.0 / 573591083, 0, 0, -433636366.0 / 683701615,
       -421739975.0 / 2616292301, 100302831.0 / 723423059,
       790204164.0 / 839813087, 800635310.0 / 3783071287, 0, 0, 0, 0, 0},
      {246121993.0 / 1340847787, 0, 0, -37695042795.0 / 15268766246,
       -309121744.0 / 1061227803, -12992083.0 / 490766935,
       6005943493.0 / 2108947869, 393006217.0 / 1396673457,
       123872331.0 / 1001029789, 0, 0, 0, 0},
      {-1028468189.0 / 846180014, 0, 0, 8478235783.0 / 508512852,
       1311729495.0 / 1432422823, -10304129995.0 / 1701304382,
       -48777925059.0 / 3047939560, 15336726248.0 / 1032824649,
       -45442868181.0 / 3398467696, 3065993473.0 / 597172653, 0, 0, 0},
      {185892177.0 / 718116043, 0, 0, -3185094517.0 / 667107341,
       -477755414.0 / 1098053517, -703635378.0 / 230739211,
       5731566787.0 / 1027545527, 5232866602.0 / 850066563,
       -4093664535.0 / 808688257, 3962137247.0 / 1805957418,
       65686358.0 / 487910083, 0, 0},
      {403863854.0 / 491063109, 0, 0, -5068492393.0 / 434740067,
       -411421997.0 / 543043805, 652783627.0 / 914296604,
       11173962825.0 / 925320556, -13158990841.0 / 6184727034,
       3936647629.0 / 1978049680, -160528059.0 / 685178525,
       248638103.0 / 1413531060, 0, 0}};

  static constexpr long double b[] = {14005451.0 / 335480064,
                               0,
                               0,
                               0,
                               0,
                               -59238493.0 / 1068277825,
                               181606767.0 / 758867731,
                               561292985.0 / 797845732,
                               -1041891430.0 / 1371343529,
                               760417239.0 / 1151165299,
                               118820643.0 / 751138087,
                               -528747749.0 / 2220607170,
                               1.0 / 4};

  static constexpr long double c[] = {0,
                               1.0 / 18,
                               1.0 / 12,
                               1.0 / 8,
                               5.0 / 16,
                               3.0 / 8,
                               59.0 / 400,
                               93.0 / 200,
                               5490023248.0 / 9719169821,
                               13.0 / 20,
                               1201146811.0 / 1299019798,
                               1,
                               1};
};

static constexpr long double sqrt2 = 1.41421356237309514547462185873883e+00;

template <>
struct ButcherTableau<RK::DIRK, RK::SDIRK2> {
  static constexpr int order = 2;
  static constexpr int nStages = 2;
  static constexpr long double a[][2] = {{1 - sqrt2 / 2, 0},
                                  {sqrt2 / 2.0, 1 - sqrt2 / 2}};
  static constexpr long double b[] = {sqrt2 / 2, 1 - sqrt2 / 2};
  static constexpr long double c[] = {1 - sqrt2 / 2, 1};
};

template <>
struct ButcherTableau<RK::DIRK, RK::ESDIRK4> {
  static constexpr int order = 4;
  static constexpr int nStages = 6;
  static constexpr long double a[][6] = {
      {0, 0, 0, 0, 0, 0},
      {1.0 / 4, 1.0 / 4, 0, 0, 0, 0},
      {(1 - sqrt2) / 8, (1 - sqrt2) / 8, 1.0 / 4, 0, 0, 0},
      {(5 - 7 * sqrt2) / 64, (5 - 7 * sqrt2) / 64, 7 * (1 + sqrt2) / 32,
       1.0 / 4, 0, 0},
      {(-13796 - 54539 * sqrt2) / 125000, (-13796 - 54539 * sqrt2) / 125000,
       (506605 + 132109 * sqrt2) / 437500, 166 * (-97 + 376 * sqrt2) / 109375,
       1.0 / 4, 0},
      {(1181 - 987 * sqrt2) / 13782, (1181 - 987 * sqrt2) / 13782,
       47 * (-267 + 1783 * sqrt2) / 273343,
       -16 * (-22922 + 3525 * sqrt2) / 571953,
       -15625 * (97 + 376 * sqrt2) / 90749876, 1.0 / 4}};
  static constexpr long double b[] = {(1181 - 987 * sqrt2) / 13782,
                               (1181 - 987 * sqrt2) / 13782,
                               47 * (-267 + 1783 * sqrt2) / 273343,
                               -16 * (-22922 + 3525 * sqrt2) / 571953,
                               -15625 * (97 + 376 * sqrt2) / 90749876,
                               1.0 / 4};
  static constexpr long double c[] = {0,       1.0 / 2,   (2 - sqrt2) / 4,
                               5.0 / 8, 26.0 / 25, 1};
};

#endif