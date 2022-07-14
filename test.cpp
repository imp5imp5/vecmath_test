#define __forceinline inline

#include <vecmath/dag_vecMath.h>
#include <string.h>
#include <memory.h>
#include <math.h>
#include <chrono>
#include <stdio.h>
#include <unistd.h>

#ifndef G_ASSERT
#include <assert.h>
#define G_ASSERT assert
#endif

namespace mathtest
{

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
  G_ASSERT(memcmp(&a, &fv2, sizeof(a)) == 0);
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


void profile(int n)
{
  static volatile float vs = 0.0;
  static volatile vec4f vv;
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


  PROFILE_BEGIN("add")
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

  PROFILE_BEGIN("mul")
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

  PROFILE_BEGIN("rcp")
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

  PROFILE_BEGIN("rcp_x")
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

  PROFILE_BEGIN("rcp_est")
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

  PROFILE_BEGIN("rcp_est_x")
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

  PROFILE_BEGIN("div")
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

  PROFILE_BEGIN("div_x")
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
}



void base_test()
{
  TestIntVec iv;
  TestFloatVec fv;
  TestFloatVec fv2;
  vec4f a, b;
  float pos_inf = exp(1e20f);
  float neg_inf = -pos_inf;
  float nan = cos(pos_inf);

  G_ASSERT(sizeof(a) == sizeof(iv));
  G_ASSERT(sizeof(a) == sizeof(fv));

  a = v_zero();
  iv.set(0, 0, 0, 0);
  G_ASSERT(memcmp(&a, &iv, sizeof(a)) == 0);

  a = v_msbit();
  iv.set(0x80000000u, 0x80000000u, 0x80000000u, 0x80000000u);
  G_ASSERT(memcmp(&a, &iv, sizeof(a)) == 0);

  fv.set(1.0f, 2.0f, 3.0f, 4.0f);
  a = v_ldu(&fv.x);
  G_ASSERT(memcmp(&a, &fv, sizeof(a)) == 0);

  fv.set(1.0f, 2.0f, 3.0f, 4.0f);
  a = v_ldu(&fv.x);
  a = v_add(a, a);
  fv2.set(2.0f, 4.0f, 6.0f, 8.0f);
  v_stu(&fv, a);
  G_ASSERT(memcmp(&fv, &fv2, sizeof(fv)) == 0);

  fv.set(1.0f, 2.0f, 3.0f, 4.0f);
  a = v_ldu(&fv.x);
  a = v_mul(a, a);
  fv2.set(1.0f, 4.0f, 9.0f, 16.0f);
  v_stu(&fv, a);
  G_ASSERT(memcmp(&fv, &fv2, sizeof(fv)) == 0);

  fv.set(1.0f, 2.0f, 3.0f, 4.0f);
  a = v_ldu(&fv.x);
  a = v_div(a, a);
  fv2.set(1.0f, 1.0f, 1.0f, 1.0f);
  v_stu(&fv, a);
  G_ASSERT(memcmp(&fv, &fv2, sizeof(fv)) == 0);

  fv.set(-1.0f, -0.0f, -1e-30f, -1e30f);
  fv2.set(1.0f, 0.0f, 1e-30f, 1e30f);
  a = v_abs(v_ldu(&fv.x));
  G_ASSERT(memcmp(&a, &fv2, sizeof(a)) == 0);

  fv.set(neg_inf, neg_inf, neg_inf, pos_inf);
  fv2.set(pos_inf, pos_inf, pos_inf, pos_inf);
  a = v_abs(v_ldu(&fv.x));
  G_ASSERT(memcmp(&a, &fv2, sizeof(a)) == 0);

  fv.set(neg_inf, neg_inf, neg_inf, pos_inf);
  fv2.set(1e30f, 1e-30f, -1e30f, -1e-30f);
  a = v_add(v_ldu(&fv.x), v_ldu(&fv2.x));
  G_ASSERT(memcmp(&a, &fv, sizeof(a)) == 0);

  fv.set(neg_inf, neg_inf, neg_inf, pos_inf);
  fv2.set(1e30f, 1e-30f, -1e30f, -1e-30f);
  a = v_sub(v_ldu(&fv.x), v_ldu(&fv2.x));
  G_ASSERT(memcmp(&a, &fv, sizeof(a)) == 0);

  fv.set(neg_inf, pos_inf, neg_inf, pos_inf);
  fv2.set(1e30f, 1e-30f, 1e30f, 1e-30f);
  a = v_mul(v_ldu(&fv.x), v_ldu(&fv2.x));
  G_ASSERT(memcmp(&a, &fv, sizeof(a)) == 0);

  fv.set(neg_inf, neg_inf, neg_inf, pos_inf);
  fv2.set(0.0f, -0.0f, 1e30f, 0.0f);
  a = v_mul(v_ldu(&fv.x), v_ldu(&fv2.x));
  fv.set(nan, nan, neg_inf, nan);
  G_ASSERT(memcmp(&a, &fv, sizeof(a)) == 0);

  fv.set(1.25f, -2.1f, 0.f, -234.0f);
  fv2.set(1, -3, 0, -234);
  a = v_floor(v_ldu(&fv.x));
  G_ASSERT(memcmp(&a, &fv2, sizeof(a)) == 0);

  fv.set(-0.0f, -1e-30f, 1e-30f, -(1 << 23));
  fv2.set(0, -1, 0, -(1 << 23));
  a = v_floor(v_ldu(&fv.x));
  G_ASSERT(memcmp(&a, &fv2, sizeof(a)) == 0);

  fv.set(1.25f, -2.1f, 0.f, -234.0f);
  fv2.set(2, -2, 0, -234);
  a = v_ceil(v_ldu(&fv.x));
  G_ASSERT(memcmp(&a, &fv2, sizeof(a)) == 0);

  fv.set(-0.0f, -1e-6f, 1e-6f, -(1 << 23));
  fv2.set(0, 0, 1, -(1 << 23));
  a = v_ceil(v_ldu(&fv.x));
  G_ASSERT(memcmp(&a, &fv2, sizeof(a)) == 0);

  fv.set(1.0f, -2.0f, 3.0f, -4.0f);
  fv2.set(0.1f, 0.1f, 0.1f, 0.1f);
  a = v_ldu(&fv.x);
  b = v_ldu(&fv2.x);
  fv2.set(1.0f, 1.0f, 1.0f, 1.0f);
  for (int i = 0; i < 18; i++)
  {
    vec4f c = v_safediv(a, a);
    a = v_mul(a, b);
    G_ASSERT(memcmp(&c, &fv2, sizeof(a)) == 0);
  }


  fv.set(1.0f, -2.0f, 0.0f, -4.0f);
  fv2.set(0.0, -0.0, 0.0, 0.0);
  a = v_ldu(&fv.x);
  b = v_ldu(&fv2.x);
  vec4f c = v_safediv(a, b);
  G_ASSERT(isfinite( ((float*)&c)[0] ));
  G_ASSERT(isfinite( ((float*)&c)[1] ));
  G_ASSERT(isfinite( ((float*)&c)[2] ));
  G_ASSERT(isfinite( ((float*)&c)[3] ));


  vec3f pa;
  ((float*)&pa)[0] = 1;
  ((float*)&pa)[1] = 1;
  ((float*)&pa)[2] = 1;
  vec3f pb = pa;
  vec3f pp = pa;
  ((float*)&pp)[1] = 2;

  vec3f res = closest_point_in_segment(pa, pb, pp);
  G_ASSERT(A(res, 0) == 1.0f && A(res, 1) == 1.0f && A(res, 2) == 1.0f);

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
    G_ASSERT(isfinite(A(w, 0)) || A(invalid, 0));
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

    G_ASSERT(v_test_segment_box_intersection_dir(start, dir, box) == 1);
  }
}

} // namespace mathtest


int main()
{
//  mathtest::base_test();
//  mathtest::advanced_test();
  mathtest::profile(20000005);
//  mathtest::profile(1);
  return 0;
}
