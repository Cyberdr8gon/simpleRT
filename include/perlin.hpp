#pragma once

#include "vec3.hpp"
#include "random_gen.h"


inline float perlin_interp(vec3 c[2][2][2], double u, double v, double w) {
    auto uu = u*u*(3-2*u);
    auto vv = v*v*(3-2*v);
    auto ww = w*w*(3-2*w);

    auto accum = 0.0;

    for (int i=0; i < 2; i++)
        for (int j=0; j < 2; j++)
            for (int k=0; k < 2; k++) {
                vec3 weight_v(u-i, v-j, w-k);
                accum += (i*uu + (1-i)*(1-uu))*
                    (j*vv + (1-j)*(1-vv))*
                    (k*ww + (1-k)*(1-ww))*dot(c[i][j][k], weight_v);
            }

    return accum;
}


class perlin
{
public:
  float noise(const vec3& p) const {
    auto u = p.x() - floor(p.x());
    auto v = p.y() - floor(p.y());
    auto w = p.z() - floor(p.z());
    auto i = static_cast<int>(floor(p.x()));
    auto j = static_cast<int>(floor(p.y()));
    auto k = static_cast<int>(floor(p.z()));
    vec3 c[2][2][2];

    for (int di=0; di < 2; di++)
        for (int dj=0; dj < 2; dj++)
            for (int dk=0; dk < 2; dk++)
                c[di][dj][dk] = ranvec[
                    perm_x[(i+di) & 255] ^
                    perm_y[(j+dj) & 255] ^
                    perm_z[(k+dk) & 255]
                ];

    return perlin_interp(c, u, v, w); 
  }

  float turb(const vec3& p, int depth=7) const {
    float accum = 0;
    vec3 temp_p = p;
    float weight = 1.0;
    for (int i = 0; i < depth; ++i) {
      accum+= weight*noise(temp_p);
      weight *=0.5;
      temp_p *= 2;
    }
    return fabs(accum);
  }
  static vec3* ranvec;
  static int *perm_x;
  static int *perm_y;
  static int *perm_z;
};

static vec3* perlin_generate() {
  vec3* p = new vec3[256];
  for (int i = 0; i < 256; ++i) {
    p[i] = unit_vector(vec3(-1.0 + 2.0*uniform(gen), -1.0 + 2.0*uniform(gen), -1.0 + 2.0*uniform(gen)));
  }
  return p;
}


void permute(int* p, int n) {
  for (int i = n-1; i > 0; i--) {
    int target = int(uniform(gen) * (i+1));
    int tmp = p[i];
    p[i] = p[target];
    p[target] = tmp;
  }
}

static int* perline_generate_perm() {
  int* p = new int[256];
  for (int i = 0; i < 256; ++i) {
    p[i] = i;
  }
  permute(p, 256);
  return p;
}

vec3* perlin::ranvec = perlin_generate();
int* perlin::perm_x = perline_generate_perm();
int* perlin::perm_y = perline_generate_perm();
int* perlin::perm_z = perline_generate_perm();


