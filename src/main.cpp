#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <ppl.h>
#include "random_gen.h"

#define _USE_MATH_DEFINES
#include <cmath>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"


#include "vec3.hpp"
#include "ray.hpp"
#include "hitablelist.hpp"
#include "sphere.h"
#include "float.h"
#include "camera.h"
#include "bvh.hpp"
#include "texture.hpp"

void write_to_buffer(unsigned char* databuffer, size_t index, vec3 color) {
  int ir = int(255.99*color[0]);
  int ig = int(255.99*color[1]);
  int ib = int(255.99*color[2]);

  databuffer[index] = ir;
  databuffer[index+1] = ig;
  databuffer[index+2] = ib;

}

vec3 random_in_unit_sphere() {
  vec3 p;
  do {
    p = 2.0*vec3(uniform(gen), uniform(gen), uniform(gen)) - vec3(1,1,1);
  } while(p.squared_length() >= 1.0);
  return p;
}


vec3 reflect(const vec3& v, const vec3& n) {
  return v - 2*dot(v, n) * n;
}

bool refract(const vec3& v, const vec3& n, float ni_over_nt, vec3& refracted) {
  vec3 uv = unit_vector(v);
  float dt = dot(uv, n);
  float discriminant = 1.0 - ni_over_nt*ni_over_nt * (1-dt*dt);
  if(discriminant > 0) {
    refracted = ni_over_nt * (uv - n*dt) - n*sqrt(discriminant);
    return true;
  } else {
    return false;
  }

}

float schlick(float cosine, float ref_idx) {
  float r0 = (1-ref_idx) / (1+ref_idx);
  r0 = r0 * r0;
  return r0 + (1-r0)*pow((1-cosine), 5);
}

class material
{
public:
  virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const = 0;
};

class lambertian : public material {
  public:
    lambertian(texture* a) : albedo(a) {}
    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
      vec3 target = rec.p + rec.normal + random_in_unit_sphere();
      scattered = ray(rec.p, target-rec.p, r_in.time());
      attenuation = albedo->value(0,0, rec.p);
      return true;
    }
    texture* albedo;
};

class metal : public material {
public:
  metal(const vec3& a, float f) : albedo(a) {
    if (f < 1){
      fuzz = f;
    } else {
      fuzz = 1;
    }
  }
    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
        vec3 reflected = reflect(unit_vector(r_in.direction()), rec.normal);
        scattered = ray(rec.p, reflected + fuzz*random_in_unit_sphere(), r_in.time());
        attenuation = albedo;
        return (dot(scattered.direction(), rec.normal) > 0);
    }

    vec3 albedo;
    float fuzz;
};

class dielectric : public material {
  public:
    dielectric(float ri) : ref_idx(ri) {}
    virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const {
      vec3 outward_normal;
      vec3 reflected = reflect(r_in.direction(), rec.normal);
      float ni_over_nt;
      attenuation = vec3(1.0, 1.0, 1.0);
      vec3 refracted;
      float reflect_prob;
      float cosine;
      if(dot(r_in.direction(), rec.normal) > 0) {
        outward_normal = -rec.normal;
        ni_over_nt = ref_idx;
        cosine = ref_idx * dot(r_in.direction(), rec.normal) / r_in.direction().length();
      } else {
        outward_normal = rec.normal;
        ni_over_nt = 1.0 / ref_idx;
        cosine = -dot(r_in.direction(), rec.normal) / r_in.direction().length();
      }

      if(refract(r_in.direction(), outward_normal, ni_over_nt, refracted)) {
        reflect_prob = schlick(cosine, ref_idx);
      } else {
        scattered = ray(rec.p, reflected, r_in.time());
        reflect_prob = 1.0;
      }
      if(uniform(gen) < reflect_prob) {
        scattered = ray(rec.p, reflected, r_in.time());
      } else {
        scattered = ray(rec.p, refracted, r_in.time());
      }
      return true;
    }

    float ref_idx;
};






float hit_sphere(const vec3& center, float radius, const ray& r) {
  vec3 oc = r.origin() - center;
  float a = dot(r.direction(), r.direction());
  float b = 2.0 * dot(oc, r.direction());
  float c = dot(oc, oc) - radius * radius;
  float discriminant = b*b - 4*a*c;
  if(discriminant < 0) {
    return -1.0;
  } else {
    return (-b -sqrt(discriminant))/(2.0 * a);
  }
}

vec3 color(const ray& r, hitable* world, int depth) {
  hit_record rec;
  if(world->hit(r, 0.001, FLT_MAX, rec)) {
    ray scattered;
    vec3 attenuation;
    if(depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
      return attenuation*color(scattered, world, depth+1);
    } else {
      return vec3(0,0,0);
    }
  } else {
    vec3 unit_direction = unit_vector(r.direction());
    float t = 0.5*(unit_direction.y() + 1.0);
    return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
  }
}

hitable* two_spheres() {
  texture* checkerTex = new checker_texture(new constant_texture(vec3(0.2, 0.3, 0.1)), new constant_texture(vec3(0.9,0.9,0.9)));
  int n = 50;
  hitable** list = new hitable*[n+1];
  list[0] = new sphere(vec3(0, -10, 0), 10, new lambertian(checkerTex));
  list[1] = new sphere(vec3(0, 10, 0), 10, new lambertian(checkerTex));

  return new hitable_list(list, 2);
}

hitable* two_perlin_spheres() {
  texture* pertext = new noise_texture(2);
  hitable** list = new hitable*[2];
  list[0] = new sphere(vec3(0,-1000, 0), 1000, new lambertian(pertext));
  list[1] = new sphere(vec3(0,2, 0), 2, new lambertian(pertext));
  return new hitable_list(list,2);
}

hitable* random_scene(float time0, float time1) {
  int n = 50000;
  hitable** list = new hitable*[n+1];

  texture* checkerTex = new checker_texture(new constant_texture(vec3(0.2, 0.3, 0.1)), new constant_texture(vec3(0.9,0.9,0.9)));
  list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(checkerTex));
  int i = 1;
  for (int a = -11; a < 11; a++) {
    for (int b = -11; b < 11; b++) {
      float choose_mat = uniform(gen);
      vec3 center(a+0.9*uniform(gen),0.2, b+0.9*uniform(gen));
      if((center - vec3(4,0.2, 0)).length() > 0.9) {
        if(choose_mat < 0.8) { //difuse
          list[i++] = new moving_sphere(center, center+vec3(0.0, 0.5 * uniform(gen), 0), 0.0, 1.0, 0.2, new lambertian(new constant_texture(vec3(uniform(gen) * uniform(gen), uniform(gen)*uniform(gen), uniform(gen)*uniform(gen)))));
        } else if (choose_mat < 0.9) { //metal
          list[i++] = new sphere(center, 0.2,
              new metal(vec3(0.5*(1 + uniform(gen)),0.5*(1 + uniform(gen)),0.5*(1 + uniform(gen))), 0.5 * uniform(gen)));
        } else { // dielectric
          list[i++] = new sphere(center, 0.2, new dielectric(1.5));
        }
      }
    }
  }
  list[i++] = new sphere(vec3(0, 1,0), 1.0, new dielectric(1.5));
  list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new lambertian(new constant_texture(vec3(0.4, 0.2, 0.1))));
  list[i++] = new sphere(vec3(4,1,0), 1.0, new metal(vec3(0.7, 0.6, 0.5), 0.0));

  return new bvh_node(list, i, time0, time1);
}


int main(int argc, char *argv[]) {

  // setup size of render image
  //const int nx = 2400;
  //const int ny = 1200;
  const int size_multiplier = 4;
  const int nx = 200*size_multiplier;
  const int ny = 100*size_multiplier;

  // samples for Anti-aliasing
  const int ns = 100;

  const int channels_num = 3;

  unsigned char* data = static_cast<unsigned char*>(malloc(sizeof(char)*nx * ny * channels_num));

  vec3 lower_left_corner(-2.0, -1.0, -1.0);
  vec3 horizontal(4.0, 0.0, 0.0);
  vec3 vertical(0.0, 2.0, 0.0);
  vec3 origin(0.0, 0.0, 0.0);

  float R = cos(M_PI/4.0);

  float time0 = 0.0;
  float time1 = 1.0;

  //hitable* world = random_scene(time0, time1);
  //hitable* world = two_spheres();
  hitable* world = two_perlin_spheres();

  vec3 lookfrom(5,1.5,2);
  vec3 lookat(0,0,0);
  float dist_to_focus = (lookfrom-lookat).length();
  float aperature = 0.0;

  camera cam(lookfrom, lookat, vec3(0,1,0), 90, float(nx)/float(ny), aperature, dist_to_focus, time0, time1);

//  for (int j = ny-1; j >=0; --j) {
  concurrency::parallel_for(int(0), ny-1, [&](size_t j) {
    for (int i = 0; i < nx; ++i) {
      vec3 col(0,0,0);
      for(int s = 0; s < ns; s++) {
        float u = float(i + uniform(gen)) / float(nx);
        float v = float(j + uniform(gen)) / float(ny);
        ray r = cam.get_ray(u,v);
        vec3 p = r.point_at_parameter(2.0);
        col += color(r,world, 0);
      }

      // normalize summed color samples
      col /= float(ns);

      // gamma correction
      col = vec3(sqrt(col[0]), sqrt(col[1]), sqrt(col[2]));

      size_t index = ((size_t(ny) - 1 - size_t(j)) * size_t(nx) + size_t(i)) * size_t(channels_num);

      write_to_buffer(data, index, col);

    }
  });

  stbi_write_bmp("test.bmp", nx, ny, channels_num, data);
  return 0;
}
