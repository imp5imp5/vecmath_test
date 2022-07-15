#define __forceinline inline

#include <vecmath/dag_vecMath.h>
#include <string.h>
#include <memory.h>
#include <math.h>
#include <chrono>
#include <stdio.h>
#include <stdint.h>

#if defined(DAGOR_DBGLEVEL)
#define DAGOR_MATH 1
#else
#define DAGOR_MATH 0
#endif


#if DAGOR_MATH == 0

#define INLINE inline

template <typename T> INLINE T max(T a, T b)
{
  return a > b ? a : b;
}

template <typename T> INLINE T min(T a, T b)
{
  return a < b ? a : b;
}

template <typename T> INLINE T clamp(T t, const T min_val, const T max_val)
{
  DisablePointersInMath<T>();
  if (t < min_val)
    t = min_val;
  if (t <= max_val)
    return t;
  return max_val;
}

// fast versions, |Relative Error| <= 1.5 * 2^(-12)
INLINE float fastinvsqrt(float x) { return v_extract_x(v_rsqrt_fast_x(v_ldu_x(&x))); }
INLINE float fastinv(float x)     { return v_extract_x(v_rcp_est_x(v_ldu_x(&x))); }
INLINE float fastsqrt(float x)    { vec4f vx = v_ldu_x(&x); return v_extract_x(v_mul_x(vx, v_rsqrt_fast_x(vx))); }

// precise version
INLINE float invsqrt(float x)     { return v_extract_x(v_rsqrt_x(v_ldu_x(&x))); }

INLINE void sincos(float rad, float &s, float &c)
{
#if defined(_GNU_SOURCE)
  sincosf(rad, &s, &c);
#else
  s = sin(rad);
  c = cos(rad);
#endif
}


#endif




#define C_MINUS        0x80000000u

#define P_ZERO         0x00000000u
#define P_INF          0x7f800000u
#define P_NAN          0x7fc00000u
#define P_NAN_2        0x7f800001u
#define P_NAN_3        0xff800002u
#define P_NAN_4        0x7fffffffu
#define P_ONE          0x3f800000u
#define P_SMALL        0x21160f25u  // 5e-19
#define P_BIG          0x5dda5e01u  // 1.96e+18
#define P_FLT_MIN      0x00800000u
#define P_HALF_P_EPS   0x3f000001u  // 0.500000059605
#define P_HALF_M_EPS   0x3effffffu  // 0.499999970198
#define P_10           0x41200000u
#define P_10_75        0x412c0000u  // 10.75
#define P_30_125       0x41f10000u  // 30.125
#define P_BIG_INT      0x4effffffu  // 2147483520.0
#define P_DENORM       0x00000008u  // 1.12e-44
#define N_ZERO         0x80000000u
#define N_INF          (C_MINUS | 0x7f800000u)
#define N_NAN          (C_MINUS | 0x7fc00000u)
#define N_NAN_2        (C_MINUS | 0x7f800001u)
#define N_NAN_3        (C_MINUS | 0xff800001u)
#define N_NAN_4        (C_MINUS | 0x7fffffffu)
#define N_ONE          (C_MINUS | 0x3f800000u)
#define N_SMALL        (C_MINUS | 0x21160f25u)  // -5e-19
#define N_BIG          (C_MINUS | 0x5dda5e01u)  // -1.96e+18
#define N_FLT_MIN      (C_MINUS | 0x00800000u)
#define N_HALF_P_EPS   (C_MINUS | 0x3f000001u)  // -0.500000059605
#define N_HALF_M_EPS   (C_MINUS | 0x3effffffu)  // -0.499999970198
#define N_10           (C_MINUS | 0x41200000u)
#define N_10_75        (C_MINUS | 0x412c0000u)  // -10.75
#define N_30_125       (C_MINUS | 0x41f10000u)  // -30.125
#define N_BIG_INT      (C_MINUS | 0x4effffffu)  // -2147483520.0
#define N_DENORM       (C_MINUS | 0x00000008u)  // -1.12e-44

#define FINITE     0
#define INF        1
#define NAN        2
#define ANY        4
#define ANY_FINITE 5

static bool allow_denorm_values = true;

typedef float (*SingleArgFunctionPtr)(float x);


uint32_t ci(float f)
{
  uint32_t x;
  memcpy(&x, &f, sizeof(float));
  return x;
}

float cf(uint32_t i)
{
  float x;
  memcpy(&x, &i, sizeof(float));
  return x;
}


uint32_t test_inf_nan(float x)
{
  uint32_t res = 0;
  uint32_t i = ci(x);
  if ((i & 0x7f800000u) == 0x7f800000u && (i & 0x007fffffu) == 0u)
    res = INF;
  if ((i & 0x7f800000u) == 0x7f800000u && (i & 0x007fffffu) != 0u)
    res = NAN;
  return res;
}

bool is_denorm_input(uint32_t x)
{
  return (x & 0x7FFFFFFFu) <= 0x007fffffu;
}


/*
  P_ZERO       ,
  P_INF        ,
  P_NAN        ,
  P_NAN_2      ,
  P_NAN_3      ,
  P_NAN_4      ,
  P_ONE        ,
  P_SMALL      ,
  P_BIG        ,
  P_FLT_MIN    ,
  P_HALF_P_EPS ,
  P_HALF_M_EPS ,
  P_10         ,
  P_10_75      ,
  P_30_125     ,
  P_BIG_INT    ,
  P_DENORM     ,
  N_ZERO       ,
  N_INF        ,
  N_NAN        ,
  N_NAN_2      ,
  N_NAN_3      ,
  N_NAN_4      ,
  N_ONE        ,
  N_SMALL      ,
  N_BIG        ,
  N_FLT_MIN    ,
  N_HALF_P_EPS ,
  N_HALF_M_EPS ,
  N_10         ,
  N_10_75      ,
  N_30_125     ,
  N_BIG_INT    ,
  N_DENORM     ,
*/


void on_fail_test(const char * test_str)
{
  printf("FAILED: %s\n", test_str);
  // TODO: error = true;
}

bool float_rel_eq(float a, float b, float rel_threshold)
{
  if (fabs(a - b) < 1e-18f)
    return true;

  float mx = max(fabs(a), fabs(b));
  float mn = min(fabs(a), fabs(b));

  float t = mx * rel_threshold;
  return (fabs(a - b) < t);
}

#define TEST(x, str) \
  { \
    if (!(x)) \
    { \
      printf("\nat file: %s, line: %d\n", __FILE__, __LINE__); \
      on_fail_test("assertion failed: " str ", " #x); \
    } \
  }

#define TEST_EQ(func, x, expect) \
  { \
    uint32_t v = ci(func(cf(x))); \
    bool ok = (v == (expect)); \
    if ((expect == P_ZERO || expect == N_ZERO) && (v == P_ZERO || v == N_ZERO)) \
      ok = true; \
    if (!ok) \
    { \
      printf("\n%s(%g) = %g,  expected %g\n", #func, cf(x), func(cf(x)), cf(expect)); \
      on_fail_test(#func "(" #x ") != " #expect); \
    } \
  }


#define TEST_THR(func, x, expect, threshold) \
  { \
    float v = func(cf(x)); \
    bool ok = float_rel_eq(v, expect, threshold); \
    if (!ok) \
    { \
      printf("\n%s(%g) = %g,  expected %g\n", #func, cf(x), v, expect); \
      on_fail_test(#func "(" #x ") != " #expect); \
    } \
  }


#define TEST_FINITE(func, x) \
  { \
    float v = func(cf(x)); \
    bool ok = test_inf_nan(v) != INF; \
    if (!ok) \
    { \
      printf("\n%s(%g) = %g\n", #func, cf(x), v); \
      on_fail_test( #func "(" #x ") is not a finite number"); \
    } \
  }


#define TEST_INF(func, x) \
  { \
    float v = func(cf(x)); \
    bool ok = test_inf_nan(v) == INF; \
    if (!ok) \
    { \
      printf("\n%s(%g) = %g\n", #func, cf(x), v); \
      on_fail_test( #func "(" #x ") is not a finite number"); \
    } \
  }


#define TEST_NAN(func, x) \
  { \
    float v = func(cf(x)); \
    bool ok = test_inf_nan(v) == NAN; \
    if (!ok) \
    { \
      printf("\n%s(%g) = %g\n", #func, cf(x), v); \
      on_fail_test( #func "(" #x ") is not a NaN"); \
    } \
  }


void test_ceil()
{
  TEST_EQ( ceilf, P_ZERO       , P_ZERO )
  TEST_EQ( ceilf, P_ONE        , P_ONE )
  TEST_EQ( ceilf, P_SMALL      , P_ONE )
  TEST_EQ( ceilf, P_FLT_MIN    , P_ONE )
  TEST_EQ( ceilf, P_HALF_P_EPS , P_ONE )
  TEST_EQ( ceilf, P_HALF_M_EPS , P_ONE )
  TEST_EQ( ceilf, P_10         , P_10 )
  TEST_EQ( ceilf, P_10_75      , ci(11.0f) )
  TEST_EQ( ceilf, P_30_125     , ci(31.0f) )
  TEST_EQ( ceilf, P_BIG_INT    , P_BIG_INT )
  TEST_EQ( ceilf, P_DENORM     , P_ONE )
  TEST_EQ( ceilf, N_ZERO       , P_ZERO )
  TEST_EQ( ceilf, N_ONE        , N_ONE )
  TEST_EQ( ceilf, N_SMALL      , P_ZERO )
  TEST_EQ( ceilf, N_FLT_MIN    , P_ZERO )
  TEST_EQ( ceilf, N_HALF_P_EPS , P_ZERO )
  TEST_EQ( ceilf, N_HALF_M_EPS , P_ZERO )
  TEST_EQ( ceilf, N_10         , N_10 )
  TEST_EQ( ceilf, N_10_75      , N_10 )
  TEST_EQ( ceilf, N_30_125     , ci(-30.0f) )
  TEST_EQ( ceilf, N_BIG_INT    , N_BIG_INT )
  TEST_EQ( ceilf, N_DENORM     , P_ZERO )

  TEST_FINITE( ceilf, P_BIG )
  TEST_FINITE( ceilf, N_BIG )
}

void test_floor()
{
  TEST_EQ( floorf, P_ZERO       , P_ZERO )
  TEST_EQ( floorf, P_ONE        , P_ONE )
  TEST_EQ( floorf, P_SMALL      , P_ZERO )
  TEST_EQ( floorf, P_FLT_MIN    , P_ZERO )
  TEST_EQ( floorf, P_HALF_P_EPS , P_ZERO )
  TEST_EQ( floorf, P_HALF_M_EPS , P_ZERO )
  TEST_EQ( floorf, P_10         , P_10 )
  TEST_EQ( floorf, P_10_75      , P_10 )
  TEST_EQ( floorf, P_30_125     , ci(30.0f) )
  TEST_EQ( floorf, P_BIG_INT    , P_BIG_INT )
  TEST_EQ( floorf, P_DENORM     , P_ZERO )
  TEST_EQ( floorf, N_ZERO       , P_ZERO )
  TEST_EQ( floorf, N_ONE        , N_ONE )
  TEST_EQ( floorf, N_SMALL      , N_ONE )
  TEST_EQ( floorf, N_FLT_MIN    , N_ONE )
  TEST_EQ( floorf, N_HALF_P_EPS , N_ONE )
  TEST_EQ( floorf, N_HALF_M_EPS , N_ONE )
  TEST_EQ( floorf, N_10         , N_10 )
  TEST_EQ( floorf, N_10_75      , ci(-11.0f) )
  TEST_EQ( floorf, N_30_125     , ci(-31.0f) )
  TEST_EQ( floorf, N_BIG_INT    , N_BIG_INT )
  TEST_EQ( floorf, N_DENORM     , N_ONE )

  TEST_FINITE( floorf, P_BIG )
  TEST_FINITE( floorf, N_BIG )
}

void test_fastinvsqrt()
{
  TEST_EQ( fastinvsqrt, P_ZERO       , P_INF )
  TEST_EQ( fastinvsqrt, P_INF        , P_ZERO )

  TEST_FINITE( fastinvsqrt, P_SMALL )
  TEST_FINITE( fastinvsqrt, P_FLT_MIN )

  TEST_THR( fastinvsqrt, P_ONE        , 1.0 / sqrt(cf(P_ONE)), 0.001f )
  TEST_THR( fastinvsqrt, P_BIG        , 1.0 / sqrt(cf(P_BIG)), 0.001f )
  TEST_THR( fastinvsqrt, P_HALF_P_EPS , 1.0 / sqrt(cf(P_HALF_P_EPS)), 0.001f )
  TEST_THR( fastinvsqrt, P_HALF_M_EPS , 1.0 / sqrt(cf(P_HALF_M_EPS)), 0.001f )
  TEST_THR( fastinvsqrt, P_10         , 1.0 / sqrt(cf(P_10)), 0.001f )
  TEST_THR( fastinvsqrt, P_10_75      , 1.0 / sqrt(cf(P_10_75)), 0.001f )
  TEST_THR( fastinvsqrt, P_30_125     , 1.0 / sqrt(cf(P_30_125)), 0.001f )
  TEST_THR( fastinvsqrt, P_BIG_INT    , 1.0 / sqrt(cf(P_BIG_INT)), 0.001f )
}

void test_invsqrt()
{
  TEST_EQ( invsqrt, P_ZERO       , P_INF )
  TEST_EQ( invsqrt, P_INF        , P_ZERO )

  TEST_FINITE( invsqrt, P_SMALL )
  TEST_FINITE( invsqrt, P_FLT_MIN )

  TEST_THR( invsqrt, P_ONE        , 1.0 / sqrt(cf(P_ONE)), 0.00001f )
  TEST_THR( invsqrt, P_BIG        , 1.0 / sqrt(cf(P_BIG)), 0.00001f )
  TEST_THR( invsqrt, P_HALF_P_EPS , 1.0 / sqrt(cf(P_HALF_P_EPS)), 0.00001f )
  TEST_THR( invsqrt, P_HALF_M_EPS , 1.0 / sqrt(cf(P_HALF_M_EPS)), 0.00001f )
  TEST_THR( invsqrt, P_10         , 1.0 / sqrt(cf(P_10)), 0.00001f )
  TEST_THR( invsqrt, P_10_75      , 1.0 / sqrt(cf(P_10_75)), 0.00001f )
  TEST_THR( invsqrt, P_30_125     , 1.0 / sqrt(cf(P_30_125)), 0.00001f )
  TEST_THR( invsqrt, P_BIG_INT    , 1.0 / sqrt(cf(P_BIG_INT)), 0.00001f )
}

void test_fastsqrt()
{
  TEST_THR( fastsqrt, P_ZERO       , sqrt(cf( P_ZERO       )), 0.001f )
  TEST_THR( fastsqrt, P_ONE        , sqrt(cf( P_ONE        )), 0.001f )
  TEST_THR( fastsqrt, P_SMALL      , sqrt(cf( P_SMALL      )), 0.001f )
  TEST_THR( fastsqrt, P_BIG        , sqrt(cf( P_BIG        )), 0.001f )
  TEST_THR( fastsqrt, P_FLT_MIN    , sqrt(cf( P_FLT_MIN    )), 0.001f )
  TEST_THR( fastsqrt, P_HALF_P_EPS , sqrt(cf( P_HALF_P_EPS )), 0.001f )
  TEST_THR( fastsqrt, P_HALF_M_EPS , sqrt(cf( P_HALF_M_EPS )), 0.001f )
  TEST_THR( fastsqrt, P_10         , sqrt(cf( P_10         )), 0.001f )
  TEST_THR( fastsqrt, P_10_75      , sqrt(cf( P_10_75      )), 0.001f )
  TEST_THR( fastsqrt, P_30_125     , sqrt(cf( P_30_125     )), 0.001f )
  TEST_THR( fastsqrt, P_BIG_INT    , sqrt(cf( P_BIG_INT    )), 0.001f )
  TEST_THR( fastsqrt, P_DENORM     , sqrt(cf( P_DENORM     )), 0.001f )
  TEST_THR( fastsqrt, N_ZERO       , sqrt(cf( N_ZERO       )), 0.001f )
}

void test_fastinv()
{
  TEST_INF( fastinv, P_ZERO   )
  TEST_INF( fastinv, N_ZERO   )
  TEST_INF( fastinv, P_DENORM )
  TEST_INF( fastinv, N_DENORM )

  TEST_THR( fastinv, P_INF        , 1.0 / (cf( P_INF        )), 0.001f )
  TEST_THR( fastinv, P_ONE        , 1.0 / (cf( P_ONE        )), 0.001f )
  TEST_THR( fastinv, P_SMALL      , 1.0 / (cf( P_SMALL      )), 0.001f )
  TEST_THR( fastinv, P_BIG        , 1.0 / (cf( P_BIG        )), 0.001f )
  TEST_THR( fastinv, P_FLT_MIN    , 1.0 / (cf( P_FLT_MIN    )), 0.001f )
  TEST_THR( fastinv, P_HALF_P_EPS , 1.0 / (cf( P_HALF_P_EPS )), 0.001f )
  TEST_THR( fastinv, P_HALF_M_EPS , 1.0 / (cf( P_HALF_M_EPS )), 0.001f )
  TEST_THR( fastinv, P_10         , 1.0 / (cf( P_10         )), 0.001f )
  TEST_THR( fastinv, P_10_75      , 1.0 / (cf( P_10_75      )), 0.001f )
  TEST_THR( fastinv, P_30_125     , 1.0 / (cf( P_30_125     )), 0.001f )
  TEST_THR( fastinv, P_BIG_INT    , 1.0 / (cf( P_BIG_INT    )), 0.001f )
  TEST_THR( fastinv, N_INF        , 1.0 / (cf( N_INF        )), 0.001f )
  TEST_THR( fastinv, N_ONE        , 1.0 / (cf( N_ONE        )), 0.001f )
  TEST_THR( fastinv, N_SMALL      , 1.0 / (cf( N_SMALL      )), 0.001f )
  TEST_THR( fastinv, N_BIG        , 1.0 / (cf( N_BIG        )), 0.001f )
  TEST_THR( fastinv, N_FLT_MIN    , 1.0 / (cf( N_FLT_MIN    )), 0.001f )
  TEST_THR( fastinv, N_HALF_P_EPS , 1.0 / (cf( N_HALF_P_EPS )), 0.001f )
  TEST_THR( fastinv, N_HALF_M_EPS , 1.0 / (cf( N_HALF_M_EPS )), 0.001f )
  TEST_THR( fastinv, N_10         , 1.0 / (cf( N_10         )), 0.001f )
  TEST_THR( fastinv, N_10_75      , 1.0 / (cf( N_10_75      )), 0.001f )
  TEST_THR( fastinv, N_30_125     , 1.0 / (cf( N_30_125     )), 0.001f )
  TEST_THR( fastinv, N_BIG_INT    , 1.0 / (cf( N_BIG_INT    )), 0.001f )
}


void test_sincos()
{
  for (int i = -10; i <= 10; i++)
  {
    for (double f = -0.001; f < 0.001; f += 0.0000001)
    {
      float a = float(3.14159265358979323846 * 0.5 * i + f);
      float s;
      float c;
      sincos(a, s, c);

      if (fabs(s) > 1.0f)
      {
        printf("\nsincos(%0.9g), sin = %0.9g\n", a, s);
        on_fail_test("sincos > 1");
      }

      if (fabs(c) > 1.0f)
      {
        printf("\nsincos(%0.9g), cos = %0.9g\n", a, c);
        on_fail_test("sincos > 1");
      }
    }
  }

  float arg = -1.0f;
  for (int i = 0; i < 20; i++, arg *= 10)
  {
    float s;
    float c;
    sincos(arg, s, c);
    bool ok = (s >= -1.0f && s <= 1.0f && c >= -1.0f && c <= 1.0f);
    if (!ok)
    {
      printf("\nsincos(%0.9g), sin = %0.9g, cos = %0.9g\n", arg, s, c);
      on_fail_test("sincos out of range");
    }
  }

}


template <int N>
float removeNan(float x) { return v_extract(v_remove_nan(v_splats(x)), N); }

void test_v_remove_nan()
{
  #define TEST_REMOVE_NAN(i) \
    TEST_EQ( removeNan<i>, P_ZERO       , P_ZERO ) \
    TEST_EQ( removeNan<i>, P_INF        , P_INF ) \
    TEST_EQ( removeNan<i>, P_NAN        , P_ZERO ) \
    TEST_EQ( removeNan<i>, P_NAN_2      , P_ZERO ) \
    TEST_EQ( removeNan<i>, P_NAN_3      , P_ZERO ) \
    TEST_EQ( removeNan<i>, P_NAN_4      , P_ZERO ) \
    TEST_EQ( removeNan<i>, P_ONE        , P_ONE ) \
    TEST_EQ( removeNan<i>, P_FLT_MIN    , P_FLT_MIN ) \
    TEST_EQ( removeNan<i>, P_HALF_P_EPS , P_HALF_P_EPS ) \
    TEST_EQ( removeNan<i>, P_DENORM     , P_DENORM ) \
    TEST_EQ( removeNan<i>, N_ZERO       , N_ZERO ) \
    TEST_EQ( removeNan<i>, N_INF        , N_INF ) \
    TEST_EQ( removeNan<i>, N_NAN        , P_ZERO ) \
    TEST_EQ( removeNan<i>, N_NAN_2      , P_ZERO ) \
    TEST_EQ( removeNan<i>, N_NAN_3      , P_ZERO ) \
    TEST_EQ( removeNan<i>, N_NAN_4      , P_ZERO ) \
    TEST_EQ( removeNan<i>, N_ONE        , N_ONE ) \
    TEST_EQ( removeNan<i>, N_FLT_MIN    , N_FLT_MIN ) \
    TEST_EQ( removeNan<i>, N_HALF_M_EPS , N_HALF_M_EPS ) \
    TEST_EQ( removeNan<i>, N_DENORM     , N_DENORM )

  TEST_REMOVE_NAN(0)
  TEST_REMOVE_NAN(1)
  TEST_REMOVE_NAN(2)
  TEST_REMOVE_NAN(3)

  #undef TEST_REMOVE_NAN
}


//VECMATH_FINLINE vec4f VECTORCALL v_remove_not_finite_2(vec4f a)
//{
//  static vec4i_const infMask = { 0x7F800000, 0x7F800000, 0x7F800000, 0x7F800000 };
//  vec4i ai = v_cast_vec4i(a);
//  return v_cast_vec4f(v_andnoti(v_cmp_eqi(v_andi(ai, infMask), infMask), ai));
//}

template <int N>
float removeNotFinite(float x) { return v_extract(v_remove_not_finite(v_splats(x)), N); }

void test_v_remove_not_finite()
{
  #define TEST_REMOVE_NOT_FINITE(i) \
    TEST_EQ( removeNotFinite<i>, P_ZERO       , P_ZERO ) \
    TEST_EQ( removeNotFinite<i>, P_INF        , P_ZERO ) \
    TEST_EQ( removeNotFinite<i>, P_NAN        , P_ZERO ) \
    TEST_EQ( removeNotFinite<i>, P_NAN_2      , P_ZERO ) \
    TEST_EQ( removeNotFinite<i>, P_NAN_3      , P_ZERO ) \
    TEST_EQ( removeNotFinite<i>, P_NAN_4      , P_ZERO ) \
    TEST_EQ( removeNotFinite<i>, P_ONE        , P_ONE ) \
    TEST_EQ( removeNotFinite<i>, P_FLT_MIN    , P_FLT_MIN ) \
    TEST_EQ( removeNotFinite<i>, P_HALF_P_EPS , P_HALF_P_EPS ) \
    TEST_EQ( removeNotFinite<i>, P_DENORM     , P_DENORM ) \
    TEST_EQ( removeNotFinite<i>, N_ZERO       , N_ZERO ) \
    TEST_EQ( removeNotFinite<i>, N_INF        , P_ZERO ) \
    TEST_EQ( removeNotFinite<i>, N_NAN        , P_ZERO ) \
    TEST_EQ( removeNotFinite<i>, N_NAN_2      , P_ZERO ) \
    TEST_EQ( removeNotFinite<i>, N_NAN_3      , P_ZERO ) \
    TEST_EQ( removeNotFinite<i>, N_NAN_4      , P_ZERO ) \
    TEST_EQ( removeNotFinite<i>, N_ONE        , N_ONE ) \
    TEST_EQ( removeNotFinite<i>, N_FLT_MIN    , N_FLT_MIN ) \
    TEST_EQ( removeNotFinite<i>, N_HALF_M_EPS , N_HALF_M_EPS ) \
    TEST_EQ( removeNotFinite<i>, N_DENORM     , N_DENORM )

  TEST_REMOVE_NOT_FINITE(0)
  TEST_REMOVE_NOT_FINITE(1)
  TEST_REMOVE_NOT_FINITE(2)
  TEST_REMOVE_NOT_FINITE(3)

  #undef TEST_REMOVE_NOT_FINITE
}







struct TestIntVec
{
  uint32_t x, y, z, w;
  void set(uint32_t x_, uint32_t y_, uint32_t z_, uint32_t w_)
  {
    x = x_;
    y = y_;
    z = z_;
    w = w_;
  }
};

struct TestFloatVec
{
  float x, y, z, w;
  void set(float x_, float y_, float z_, float w_)
  {
    x = x_;
    y = y_;
    z = z_;
    w = w_;
  }
};

void advanced_test()
{
  TestIntVec iv;
  TestFloatVec fv;
  TestFloatVec fv2;
  vec4f a, b;
  float pos_inf = exp(1e20f);
  float neg_inf = -pos_inf;
  float nan = cos(pos_inf);

  fv.set(-0.0f, -1e-30f, 1e-30f, (1 << 23));
  fv2.set(0, 0, 1, (1 << 23));
  a = v_ceil(v_ldu(&fv.x));
  TEST(memcmp(&a, &fv2, sizeof(a)) == 0, "v_ceil");
}



#define A(x, n)   ((float*)&x)[n]
#define SET3(x, v0, v1, v2) {A(x,0)=v0; A(x,1)=v1; A(x,2)=v2;}
#define SET4(x, v0, v1, v2, v3) {A(x,0)=v0; A(x,1)=v1; A(x,2)=v2; A(x,3)=v3;}

#define PROFILE_BEGIN(name_) \
{ \
  auto t0 = std::chrono::steady_clock::now(); \
  volatile char name[] = name_;


#define PROFILE_END() \
  auto t1 = std::chrono::steady_clock::now(); \
  int tt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count(); \
  if (baseline == 0) \
    baseline = tt; \
  else \
    tt -= baseline; \
  printf("%s = %d ms\n", name, tt); \
}


float vs = 0.0;
vec4f vv;


void profile(int n)
{
  int baseline = 0;
  TestFloatVec tt;

  vec4f t;
  SET4(t, 1e-4, -1e-3, -1e-6, 1e-4);

  PROFILE_BEGIN("warm-up")
  float s = 0;
  for (int i = 0; i < 18000000; i++)
    s += acos(i * 1e-20f);
  vs = s;
  PROFILE_END()

  baseline = 0;

  PROFILE_BEGIN("baseline")
  vec4f a;
  vec4f b;
  SET4(a, 1.01, 2.01, -21.01, 0.53f);
  b = a;
  for (int i = 0; i < n; i++)
  {
    b = v_rot_1(b);
  }
  vv = v_add(vv, b);
  PROFILE_END()


  PROFILE_BEGIN("v_add")
  vec4f a;
  vec4f b;
  SET4(a, 1.001, -1.002, 0.999, -0.998);
  b = a;
  for (int i = 0; i < n; i++)
  {
    b = v_rot_1(b);
    b = v_add(a, b);
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, b);
  PROFILE_END()

  PROFILE_BEGIN("v_mul")
  vec4f a;
  vec4f b;
  SET4(a, 1.001, -1.002, 0.999, -0.998);
  b = a;
  for (int i = 0; i < n; i++)
  {
    b = v_rot_1(b);
    b = v_mul(a, b);
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, b);
  PROFILE_END()

  PROFILE_BEGIN("v_rcp")
  vec4f a;
  vec4f b;
  SET4(a, 1.001, -1.002, 0.999, -0.998);
  b = a;
  for (int i = 0; i < n; i++)
  {
    b = v_rot_1(b);
    b = v_rcp(b);
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, b);
  PROFILE_END()

  PROFILE_BEGIN("v_rcp_x")
  vec4f a;
  vec4f b;
  SET4(a, 1.001, -1.002, 0.999, -0.998);
  b = a;
  for (int i = 0; i < n; i++)
  {
    b = v_rot_1(b);
    b = v_rcp_x(b);
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, b);
  PROFILE_END()

  PROFILE_BEGIN("v_rcp_est")
  vec4f a;
  vec4f b;
  SET4(a, 1.001, -1.002, 0.999, -0.998);
  b = a;
  for (int i = 0; i < n; i++)
  {
    b = v_rot_1(b);
    b = v_rcp_est(b);
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, b);
  PROFILE_END()

  PROFILE_BEGIN("v_rcp_est_x")
  vec4f a;
  vec4f b;
  SET4(a, 1.001, -1.002, 0.999, -0.998);
  b = a;
  for (int i = 0; i < n; i++)
  {
    b = v_rot_1(b);
    b = v_rcp_est_x(b);
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, b);
  PROFILE_END()

  PROFILE_BEGIN("v_div")
  vec4f a;
  vec4f b;
  SET4(a, 1.00001, -1.00002, 0.99999, -0.99998);
  b = a;
  for (int i = 0; i < n; i++)
  {
    b = v_rot_1(b);
    b = v_div(b, a);
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, b);
  PROFILE_END()

  PROFILE_BEGIN("v_div_x")
  vec4f a;
  vec4f b;
  SET4(a, 1.00001, -1.00002, 0.99999, -0.99998);
  b = a;
  for (int i = 0; i < n; i++)
  {
    b = v_rot_1(b);
    b = v_div_x(b, a);
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, b);
  PROFILE_END()

  PROFILE_BEGIN("1/x")
  vec4f a;
  vec4f b;
  SET4(a, 1.00001, -1.00002, 0.99999, -0.99998);
  b = a;
  float xx = 1.56f;

  for (int i = 0; i < n; i++)
  {
    xx = 1.f / xx + xx * 1e-20f;
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, v_splats(xx));
  PROFILE_END()

  PROFILE_BEGIN("fastinvsqrt")
  vec4f a;
  vec4f b;
  SET4(a, 1.00001, -1.00002, 0.99999, -0.99998);
  b = a;
  float xx = 1.56f;

  for (int i = 0; i < n; i++)
  {
    xx = fastinvsqrt(xx);
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, v_splats(xx));
  PROFILE_END()

  PROFILE_BEGIN("fastinv")
  vec4f a;
  vec4f b;
  SET4(a, 1.00001, -1.00002, 0.99999, -0.99998);
  b = a;
  float xx = 1.56f;

  for (int i = 0; i < n; i++)
  {
    xx = fastinv(xx) + xx * 1e-20f;
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, v_splats(xx));
  PROFILE_END()

  PROFILE_BEGIN("fastsqrt")
  vec4f a;
  vec4f b;
  SET4(a, 1.00001, -1.00002, 0.99999, -0.99998);
  b = a;
  float xx = 1.56f;

  for (int i = 0; i < n; i++)
  {
    xx = fastsqrt(xx);
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, v_splats(xx));
  PROFILE_END()

  PROFILE_BEGIN("invsqrt")
  vec4f a;
  vec4f b;
  SET4(a, 1.00001, -1.00002, 0.99999, -0.99998);
  b = a;
  float xx = 1.56f;

  for (int i = 0; i < n; i++)
  {
    xx = invsqrt(xx);
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, v_splats(xx));
  PROFILE_END()

  PROFILE_BEGIN("sqrtf")
  vec4f a;
  vec4f b;
  SET4(a, 1.00001, -1.00002, 0.99999, -0.99998);
  b = a;
  float xx = 1.56f;

  for (int i = 0; i < n; i++)
  {
    xx = sqrtf(xx);
  }
 // printf("%g %g %g %g\n", ((float*)(&b))[0], ((float*)(&b))[1], ((float*)(&b))[2], ((float*)(&b))[3]);
  vv = v_add(vv, v_splats(xx));
  PROFILE_END()
}



void base_test()
{
  test_ceil();
  test_floor();
  test_invsqrt();
  test_fastinvsqrt();
  test_fastsqrt();
  test_fastinv();
  test_sincos();
  test_v_remove_nan();
  test_v_remove_not_finite();



  TestIntVec iv;
  TestFloatVec fv;
  TestFloatVec fv2;
  vec4f a, b;
  float pos_inf = exp(1e20f);
  float neg_inf = -pos_inf;
  float nan = cos(pos_inf);

  TEST(sizeof(a) == sizeof(iv));
  TEST(sizeof(a) == sizeof(fv));

  a = v_zero();
  iv.set(0, 0, 0, 0);
  TEST(memcmp(&a, &iv, sizeof(a)) == 0, "v_zero");

  a = v_msbit();
  iv.set(0x80000000u, 0x80000000u, 0x80000000u, 0x80000000u);
  TEST(memcmp(&a, &iv, sizeof(a)) == 0, "v_msbit");

  fv.set(1.0f, 2.0f, 3.0f, 4.0f);
  a = v_ldu(&fv.x);
  TEST(memcmp(&a, &fv, sizeof(a)) == 0, "v_ldu");

  fv.set(1.0f, 2.0f, 3.0f, 4.0f);
  a = v_ldu(&fv.x);
  a = v_add(a, a);
  fv2.set(2.0f, 4.0f, 6.0f, 8.0f);
  v_stu(&fv, a);
  TEST(memcmp(&fv, &fv2, sizeof(fv)) == 0, "v_add");

  fv.set(1.0f, 2.0f, 3.0f, 4.0f);
  a = v_ldu(&fv.x);
  a = v_mul(a, a);
  fv2.set(1.0f, 4.0f, 9.0f, 16.0f);
  v_stu(&fv, a);
  TEST(memcmp(&fv, &fv2, sizeof(fv)) == 0, "v_mul");

  fv.set(1.0f, 2.0f, 3.0f, 4.0f);
  a = v_ldu(&fv.x);
  a = v_div(a, a);
  fv2.set(1.0f, 1.0f, 1.0f, 1.0f);
  v_stu(&fv, a);
  TEST(memcmp(&fv, &fv2, sizeof(fv)) == 0, "v_div");

  fv.set(-1.0f, -0.0f, -1e-30f, -1e30f);
  fv2.set(1.0f, 0.0f, 1e-30f, 1e30f);
  a = v_abs(v_ldu(&fv.x));
  TEST(memcmp(&a, &fv2, sizeof(a)) == 0, "v_abs");

  fv.set(neg_inf, neg_inf, neg_inf, pos_inf);
  fv2.set(pos_inf, pos_inf, pos_inf, pos_inf);
  a = v_abs(v_ldu(&fv.x));
  TEST(memcmp(&a, &fv2, sizeof(a)) == 0, "v_abs inf");

  fv.set(neg_inf, neg_inf, neg_inf, pos_inf);
  fv2.set(1e30f, 1e-30f, -1e30f, -1e-30f);
  a = v_add(v_ldu(&fv.x), v_ldu(&fv2.x));
  TEST(memcmp(&a, &fv, sizeof(a)) == 0, "v_add inf");

  fv.set(neg_inf, neg_inf, neg_inf, pos_inf);
  fv2.set(1e30f, 1e-30f, -1e30f, -1e-30f);
  a = v_sub(v_ldu(&fv.x), v_ldu(&fv2.x));
  TEST(memcmp(&a, &fv, sizeof(a)) == 0, "v_sub");

  fv.set(neg_inf, pos_inf, neg_inf, pos_inf);
  fv2.set(1e30f, 1e-30f, 1e30f, 1e-30f);
  a = v_mul(v_ldu(&fv.x), v_ldu(&fv2.x));
  TEST(memcmp(&a, &fv, sizeof(a)) == 0, "v_mul");

  fv.set(neg_inf, neg_inf, neg_inf, pos_inf);
  fv2.set(0.0f, -0.0f, 1e30f, 0.0f);
  a = v_mul(v_ldu(&fv.x), v_ldu(&fv2.x));
  fv.set(nan, nan, neg_inf, nan);
  TEST(memcmp(&a, &fv, sizeof(a)) == 0, "v_mul");

  fv.set(1.25f, -2.1f, 0.f, -234.0f);
  fv2.set(1, -3, 0, -234);
  a = v_floor(v_ldu(&fv.x));
  TEST(memcmp(&a, &fv2, sizeof(a)) == 0, "v_floor");

  fv.set(-0.0f, -1e-30f, 1e-30f, -(1 << 23));
  fv2.set(0, -1, 0, -(1 << 23));
  a = v_floor(v_ldu(&fv.x));
  TEST(memcmp(&a, &fv2, sizeof(a)) == 0, "v_floor");

  fv.set(1.25f, -2.1f, 0.f, -234.0f);
  fv2.set(2, -2, 0, -234);
  a = v_ceil(v_ldu(&fv.x));
  TEST(memcmp(&a, &fv2, sizeof(a)) == 0, "v_ceil");

  fv.set(-0.0f, -1e-6f, 1e-6f, -(1 << 23));
  fv2.set(0, 0, 1, -(1 << 23));
  a = v_ceil(v_ldu(&fv.x));
  TEST(memcmp(&a, &fv2, sizeof(a)) == 0, "v_ceil");

  /*
  fv.set(1.0f, -2.0f, 3.0f, -4.0f);
  fv2.set(0.1f, 0.1f, 0.1f, 0.1f);
  a = v_ldu(&fv.x);
  b = v_ldu(&fv2.x);
  fv2.set(1.0f, 1.0f, 1.0f, 1.0f);
  for (int i = 0; i < 18; i++)
  {
    vec4f c = v_safediv(a, a);
    a = v_mul(a, b);
    TEST(memcmp(&c, &fv2, sizeof(a)) == 0, "v_safediv, mul");
  }*/


  fv.set(1.0f, -2.0f, 0.0f, -4.0f);
  fv2.set(0.0, -0.0, 0.0, 0.0);
  a = v_ldu(&fv.x);
  b = v_ldu(&fv2.x);
  vec4f c = v_safediv(a, b);
  TEST(isfinite( ((float*)&c)[0] ), "v_safediv x");
  TEST(isfinite( ((float*)&c)[1] ), "v_safediv y");
  TEST(isfinite( ((float*)&c)[2] ), "v_safediv z");
  TEST(isfinite( ((float*)&c)[3] ), "v_safediv w");


  vec3f pa;
  ((float*)&pa)[0] = 1;
  ((float*)&pa)[1] = 1;
  ((float*)&pa)[2] = 1;
  vec3f pb = pa;
  vec3f pp = pa;
  ((float*)&pp)[1] = 2;

  vec3f res = closest_point_in_segment(pa, pb, pp);
  TEST(A(res, 0) == 1.0f && A(res, 1) == 1.0f && A(res, 2) == 1.0f, "");

  {
    vec4f invalid;  // (0, 1, 0, 1) (0, 1, 1e-15, 2) (1e-15, 1, 0, 3)
    plane3f p0;
    A(p0, 0) = 0;
    A(p0, 1) = 1;
    A(p0, 2) = 0;
    A(p0, 3) = 1;
    plane3f p1;
    A(p1, 0) = 0;
    A(p1, 1) = 1;
    A(p1, 2) = 1e-15;
    A(p1, 3) = 2;
    plane3f p2;
    A(p2, 0) = 1e-15;
    A(p2, 1) = 1;
    A(p2, 2) = 0;
    A(p2, 3) = 3;
    vec3f w = three_plane_intersection(p0, p1, p2, invalid);
    TEST(isfinite(A(w, 0)) || A(invalid, 0), "three_plane_intersection");
  }

  {
    vec3f start;
    A(start, 0) = 0;
    A(start, 1) = -0.5f;
    A(start, 2) = 0;
    vec3f dir;
    A(dir, 0) = -1e-38;
    A(dir, 1) = 1;
    A(dir, 2) = 1e-38;

    bbox3f box;
    A(box.bmax, 0) = 1;
    A(box.bmax, 1) = 1;
    A(box.bmax, 2) = 1;

    A(box.bmin, 0) = -1;
    A(box.bmin, 1) = -1;
    A(box.bmin, 2) = -1;

    TEST(v_test_segment_box_intersection_dir(start, dir, box) == 1, "");
  }
}



int main()
{
  base_test();
//  mathtest::advanced_test();
  profile(20000005);
//  mathtest::profile(1);
  return 0;
}
