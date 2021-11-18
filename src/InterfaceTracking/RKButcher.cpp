#include "RKButcher.h"


constexpr int ButcherTableau<ERK, ForwardEuler>::nStages;
constexpr Real ButcherTableau<ERK, ForwardEuler>::a[][ButcherTableau<ERK, ForwardEuler>::nStages];
constexpr Real ButcherTableau<ERK, ForwardEuler>::b[];
constexpr Real ButcherTableau<ERK, ForwardEuler>::c[];

constexpr int ButcherTableau<ERK, ClassicRK4>::nStages;
constexpr Real ButcherTableau<ERK, ClassicRK4>::a[][ButcherTableau<ERK, ClassicRK4>::nStages];
constexpr Real ButcherTableau<ERK, ClassicRK4>::b[];
constexpr Real ButcherTableau<ERK, ClassicRK4>::c[];
