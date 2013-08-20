#ifndef _IBL_H_
#define _IBL_H_

#include "hdr.h"
#include "constant.h"
#include "ray.h"
#include "random.h"

#include <algorithm>
// #include <array>

class IBL {
 private:
  HDRImage hdr_;
 public:
  int width_, height_;

  inline float luminance(const Color &color) {
    return 0.298912f * color.x_ + 0.586611f * color.y_ + 0.114478f * color.z_;
  }
	

  IBL(const std::string &filename, const bool importance = false, const int importance_map_size = 0) {
    hdr_.load_unsafe(filename);
    width_ = hdr_.width();
    height_ = hdr_.height();

  }
	
  struct EnvSample {
    float u, v;
  };
  void hammersley(EnvSample *result, int n) {
    float p, u, v;
    int k, kk, pos;

    for (k=0, pos=0 ; k<n ; k++) {
      u = 0;
      for (p=0.5, kk=k ; kk ; p*=0.5, kk>>=1)
	if (kk & 1)                           // kk mod 2 == 1
	  u += p;

      v = (k + 0.5f) / n;

      result[pos].u = u;
      result[pos].v = v;
      pos ++;
    }
  }


  // DebevecのオリジナルHDR
  inline Color sample_sphere(const Vec &dir) {
    const float r = (1.0f / kPI) * acos(dir.z_) / sqrt(dir.x_ * dir.x_ + dir.y_ * dir.y_);

    float u = (dir.x_ * r + 1.0f) / 2.0f;
    float v = 1.0f - (dir.y_ * r + 1.0f) / 2.0f;
		
    if (u < 0.0f)
      u += 1.0f;
    if (v < 0.0f)
      v += 1.0f;

    const int x = (int)(u * width_) % width_;
    const int y = height_ - 1 - (int)(v * height_) % height_;

    return *hdr_.image_ptr(x, y);
  }

  // phi = [0, 2pi], theta = [0, pi]
  // theta = 0 が真上
  // v = 1が真上
  // u = [0, 1), v = [0, 1)が直接phi,thetaにマッピングされてるサンプリング方法
  inline Color sample(const Vec &dir) {
    const float theta = acos(dir.y_);
    float phi = acos(dir.x_ / sqrt(dir.x_ * dir.x_ + dir.z_ * dir.z_));
    if (dir.z_ < 0.0f)
      phi = 2.0f * kPI - phi;

    const float u = phi / (2.0f * kPI);
    const float v = 1.0f - theta / kPI;

    const int x = (int)(u * width_) % width_;
    const int y = (int)(v * height_) % height_;

    return *hdr_.image_ptr(x, y);
  }

};

#endif
