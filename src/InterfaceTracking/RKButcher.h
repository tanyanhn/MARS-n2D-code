#ifndef _RKBUTCHER_
#define _RKBUTCHER_

enum RK_Category1
{
    ERK = 1,
    DIRK,
    ARK
};

enum RK_Category2
{
    ForwardEuler = 1,
    ClassicRK4,
    nRK_Type
};

template <RK_Category1 Type1, RK_Category2 Type2>
struct ButcherTableau;

template <>
struct ButcherTableau<ERK, ForwardEuler>
{
    static constexpr int nStages = 1;
    static constexpr Real a[][1] =
        {{0}};
    static constexpr Real b[] =
        {1};
    static constexpr Real c[] =
        {0};
};

constexpr int ButcherTableau<ERK, ForwardEuler>::nStages;
constexpr Real ButcherTableau<ERK, ForwardEuler>::a[][ButcherTableau<ERK, ForwardEuler>::nStages];
constexpr Real ButcherTableau<ERK, ForwardEuler>::b[];
constexpr Real ButcherTableau<ERK, ForwardEuler>::c[];

template <>
struct ButcherTableau<ERK, ClassicRK4>
{
    static constexpr int nStages = 4;
    static constexpr Real a[][4] =
        {{0, 0, 0, 0},
         {0.5, 0, 0, 0},
         {0, 0.5, 0, 0},
         {0, 0, 1, 0}};
    static constexpr Real b[] =
        {1.0 / 6, 1.0 / 3, 1.0 / 3, 1.0 / 6};
    static constexpr Real c[] =
        {0.0, 0.5, 0.5, 1.0};
};

constexpr int ButcherTableau<ERK, ClassicRK4>::nStages;
constexpr Real ButcherTableau<ERK, ClassicRK4>::a[][ButcherTableau<ERK, ClassicRK4>::nStages];
constexpr Real ButcherTableau<ERK, ClassicRK4>::b[];
constexpr Real ButcherTableau<ERK, ClassicRK4>::c[];

#endif