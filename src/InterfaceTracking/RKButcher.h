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
struct RK_ButcherTableau;

template <>
struct RK_ButcherTableau<ERK, ForwardEuler>
{
    static constexpr int nS = 1;
    static constexpr int a[][2] =
        {{0, 1}};
    static constexpr int b[] =
        {1, 1};
    static constexpr int c[] =
        {0, 1};
};

constexpr int RK_ButcherTableau<ERK, ForwardEuler>::nS;
constexpr int RK_ButcherTableau<ERK, ForwardEuler>::a[][RK_ButcherTableau<ERK, ForwardEuler>::nS + 1];
constexpr int RK_ButcherTableau<ERK, ForwardEuler>::b[];
constexpr int RK_ButcherTableau<ERK, ForwardEuler>::c[];

template <>
struct RK_ButcherTableau<ERK, ClassicRK4>
{
    static constexpr int nS = 4;
    static constexpr int a[][5] =
        {{0, 0, 0, 0, 1},
         {1, 0, 0, 0, 2},
         {0, 1, 0, 0, 2},
         {0, 0, 1, 0, 1}};
    static constexpr int b[] =
        {1, 2, 2, 1, 6};
    static constexpr int c[] =
        {0, 1, 1, 2, 2};
};

constexpr int RK_ButcherTableau<ERK, ClassicRK4>::nS;
constexpr int RK_ButcherTableau<ERK, ClassicRK4>::a[][RK_ButcherTableau<ERK, ClassicRK4>::nS + 1];
constexpr int RK_ButcherTableau<ERK, ClassicRK4>::b[];
constexpr int RK_ButcherTableau<ERK, ClassicRK4>::c[];

#endif