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
#include "rect.h"
#include "box.hpp"
#include "float.h"
#include "camera.h"
#include "bvh.hpp"
#include "texture.hpp"
#include "material.hpp"
#include "constantMedium.hpp"

void write_to_buffer(unsigned char* databuffer, size_t index, vec3 color) {
  int ir = int(255.99*color[0]);
  int ig = int(255.99*color[1]);
  int ib = int(255.99*color[2]);
  if(ir > 255) ir = 255;
  if(ig > 255) ig = 255;
  if(ib > 255) ib = 255;

  databuffer[index] = ir;
  databuffer[index+1] = ig;
  databuffer[index+2] = ib;

}



// TODO I think remove this? replaced by hitable stuff;

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

/* 
 * color function
 */

vec3 color(const ray& r, hitable* world, int depth) {
  hit_record rec;
  if(world->hit(r, 0.001, FLT_MAX, rec)) {
    ray scattered;
    vec3 attenuation;
    vec3 emitted = rec.mat_ptr->emitted(rec.u, rec.v, rec.p);

    if(depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
      return emitted + attenuation*color(scattered, world, depth+1);
    } else {
      return emitted;
    }
  } else {
    //vec3 unit_direction = unit_vector(r.direction());
    //float t = 0.5*(unit_direction.y() + 1.0);
    //return (1.0-t)*vec3(1.0, 1.0, 1.0) + t*vec3(0.5, 0.7, 1.0);
    //
    return vec3(0,0,0);
  }
}


// ****************
// scene generators
// ****************

hitable* earth_sphere() {
  int nx, ny, nn;
  const std::string filename = "map.jpg";
  unsigned char* tex_data = stbi_load(filename.c_str(), &nx, &ny, &nn, 0);
  if (!tex_data) {
      throw std::runtime_error("File " + filename + " was not found");
  }
  material* mat = new lambertian(new image_texture(tex_data, nx, ny));

  hitable** list = new hitable*[5];
  list[0] = new sphere(vec3(0,2,0), 2, mat);
  return new hitable_list(list, 1);
}

hitable* simple_light() {
  texture* pertext = new noise_texture(4);
  hitable** list = new hitable*[4];
  list[0] = new sphere(vec3(0,-1000, 0), 1000, new lambertian(pertext));
  list[1] = new sphere(vec3(0,2, 0), 2, new lambertian(pertext));
  list[2] = new xy_rect(3,5,1,3,-2, new diffuse_light(new constant_texture(vec3(4,4,4))));
  list[3] = new sphere(vec3(0,7, 0), 2, new diffuse_light(new constant_texture(vec3(4,4,4))));
  return new hitable_list(list, 4);
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

hitable* cornell_box() {
  hitable** list = new hitable*[8];
  int i = 0;
  material* red = new lambertian(new constant_texture(vec3(0.65, 0.05, 0.05)));
  material* white = new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
  material* green = new lambertian(new constant_texture(vec3(0.12, 0.45, 0.15)));
  material* light = new diffuse_light(new constant_texture(vec3(15,15,15)));

  list[i++] = new flip_normals(new yz_rect(0,555,0, 555,555,green));
  list[i++] = new yz_rect(0,555,0, 555,0, red);
  list[i++] = new xz_rect(213,343,227, 332,554, light);
  list[i++] = new flip_normals(new xz_rect(0,555,0, 555,555, white));
  list[i++] = new xz_rect(0,555,0, 555,0, white);
  list[i++] = new flip_normals(new xy_rect(0,555,0, 555,555, white));
  list[i++] = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 165, 165), white), -18), vec3(130,0,65));
  list[i++] = new translate(new rotate_y(new box(vec3(0, 0, 0), vec3(165, 330, 165), white), 15), vec3(265, 0, 295));
  return new hitable_list(list, i);
}

hitable * cornell_smoke() {
  hitable** list = new hitable*[8];
  int i = 0;
  material* red = new lambertian(new constant_texture(vec3(0.65, 0.05, 0.05)));
  material* white = new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
  material* green = new lambertian(new constant_texture(vec3(0.12, 0.45, 0.15)));
  material* light = new diffuse_light(new constant_texture(vec3(7, 7, 7)));

  list[i++] = new flip_normals(new yz_rect(0, 555,0,555,555, green));
  list[i++] = new yz_rect(0, 555,0,555,0, red);
  list[i++] = new xz_rect(113, 443,127,432,554, light);
  list[i++] = new flip_normals(new xz_rect(0,555,0,555,555,white));
  list[i++] = new xz_rect(0,555,0,555,0,white);
  list[i++] = new flip_normals(new xy_rect(0,555,0,555,555,white));
  hitable* b1 = new translate(new rotate_y(new box(vec3(0,0,0), vec3(165, 165,165), white), -18), vec3(130, 0, 65));
  hitable* b2 = new translate(new rotate_y(new box(vec3(0,0,0), vec3(165, 330,165), white), 15), vec3(265, 0,295));

  list[i++] = new constant_medium(b1, 0.01, new constant_texture(vec3(1.0, 1.0, 1.0)));
  list[i++] = new constant_medium(b2, 0.01, new constant_texture(vec3(0.0, 0.0, 0.0)));
  return new hitable_list(list,i);
}

hitable* final() {
  int nb = 20;
  hitable** list = new hitable*[30];
  hitable** boxlist = new hitable*[10000];
  hitable** boxlist2 = new hitable*[10000];
  material* white = new lambertian(new constant_texture(vec3(0.73, 0.73, 0.73)));
  material* ground = new lambertian(new constant_texture(vec3(0.48, 0.83, 0.53)));

  int b = 0;
  for (int i = 0; i < nb; ++i) {
    for (int j = 0; j < nb; ++j) {
      float w = 100;
      float x0 = -1000 + i*w;
      float z0 = -1000 + j*w;
      float y0 = 0;
      float x1 = x0 + w;
      float y1 = 100.0*(uniform(gen) + 0.01);
      float z1 = z0 + w;
      boxlist[b++] = new box(vec3(x0, y0, z0), vec3(x1, y1, z1), ground);
    }
  }
  int l = 0;
  list[l++] = new bvh_node(boxlist, b, 0, 1);

  material* light = new diffuse_light(new constant_texture(vec3(7,7,7)));
  //list[l++] = new xz_rect(123, 423, 147, 412, 554, light);
  list[l++] = new sphere(vec3(123, 523, 554), 100, light);

  vec3 center(400, 400, 200);
  list[l++] = new moving_sphere(center, center+vec3(30,0,0), 0,1,50, new lambertian(new constant_texture(vec3(0.7, 0.3, 0.1))));

  list[l++] = new sphere(vec3(260, 150, 45), 50, new dielectric(1.5));
  list[l++] = new sphere(vec3(0, 150, 145), 50, new metal(vec3(0.8, 0.8, 0.9), 10.0));

  hitable* boundary = new sphere(vec3(360, 150, 145), 70, new dielectric(1.5));

  list[l++] = boundary;
  list[l++] = new constant_medium(boundary, 0.2, new constant_texture(vec3(0.2, 0.4, 0.9)));
  boundary = new sphere(vec3(0,0,0), 5000, new dielectric(1.5));
  list[l++] = new constant_medium(boundary, 0.0001, new constant_texture(vec3(1.0, 1.0, 1.0)));

  int nx, ny, nn;
  unsigned char* tex_data = stbi_load("map.jpg", &nx, &ny, &nn, 0);
  material* emat = new lambertian(new image_texture(tex_data, nx, ny));
  list[l++] = new sphere(vec3(400, 200, 400), 100, emat);

  texture* pertext = new noise_texture(0.1);
  list[l++] = new sphere(vec3(220, 280, 300), 800, new lambertian(pertext));
  int ns = 1000;
  for (int j = 0; j < ns; ++j) {
    boxlist2[j] = new sphere(vec3(165*uniform(gen), 165*uniform(gen), 165*uniform(gen)), 100, white);
  }
  list[l++] = new translate(new rotate_y(new bvh_node(boxlist2, ns, 0.0, 1.0), 15), vec3(-100, 270, 395));
  return new hitable_list(list,l);

}


/* 
 * main function
 */

int main(int argc, char *argv[]) {

  // setup size of render image
  //const int nx = 2400;
  //const int ny = 1200;
  const int size_multiplier = 1;
  const int nx = 200*size_multiplier;
  const int ny = 100*size_multiplier;

  // samples for Anti-aliasing
  const int ns = 10000;

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
  //hitable* world = two_perlin_spheres();
  //hitable* world = earth_sphere();
  //hitable* world = simple_light();
  //hitable* world = cornell_box();
  //hitable* world = cornell_smoke();
  hitable* world = final();

  //vec3 lookfrom(5,5,2);
  //vec3 lookat(0,2,0);
  //float dist_to_focus = (lookfrom-lookat).length();
  //float aperature = 0.0;
  //float vfov = 90;

  vec3 lookfrom(478, 278, -600);
  vec3 lookat(220, 280, 300);
  float dist_to_focus = 10.0;
  float aperature = 0.0;
  float vfov = 40.0;

  camera cam(lookfrom, lookat, vec3(0,1,0), vfov, float(nx)/float(ny), aperature, dist_to_focus, time0, time1);

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
