#ifndef _BVH_BASE_
#define _BVH_BASE_

#include <vector>
#include <algorithm>
#include <queue>

#include "vec.h"
#include "bbox.h"
#include "triangle.h"


class BVHBase {
 protected:
  struct BVHPrimitiveInfo {
    int primitiveNumber_;
    Vec centroid_;
    BBox bounds_;

  BVHPrimitiveInfo(const int pn, const BBox& b) :
    primitiveNumber_(pn), bounds_(b){
      centroid_ = 0.5f * b.pmin_ + 0.5f * b.pmax_;
    }
  };

  struct ComparePoints {
    int dim_;
  ComparePoints(const int d) : dim_(d){}
    bool operator()(const BVHPrimitiveInfo &a, const BVHPrimitiveInfo &b) const {
      return a.centroid_[dim_] < b.centroid_[dim_];
    }
  };

  void split(std::vector<BVHPrimitiveInfo> &buildData, const BBox &bbox, int start, int end, int *mid, int *dim) {
    const int nPrimitives = end - start;
    BBox centroidBounds;
    for (int i = start; i < end; ++i)
      centroidBounds = Union(centroidBounds, buildData[i].centroid_);
    *dim = centroidBounds.maximum_extent();
    *mid = (start + end) / 2;
    std::nth_element(&buildData[start], &buildData[*mid], &buildData[end-1] + 1, ComparePoints(*dim));
  }

  inline bool intersect_bbox(const BBox &bounds, const Ray &ray, const Vec& invDir, const int dirIsNeg[3]) {
    float tmin =  (bounds[  dirIsNeg[0]].x_ - ray.org_.x_) * invDir.x_;
    float tmax =  (bounds[1-dirIsNeg[0]].x_ - ray.org_.x_) * invDir.x_;
    const float tymin = (bounds[  dirIsNeg[1]].y_ - ray.org_.y_) * invDir.y_;
    const float tymax = (bounds[1-dirIsNeg[1]].y_ - ray.org_.y_) * invDir.y_;
    if ((tmin > tymax) || (tymin > tmax))
      return false;
    if (tymin > tmin) tmin = tymin;
    if (tymax < tmax) tmax = tymax;

    const float tzmin = (bounds[  dirIsNeg[2]].z_ - ray.org_.z_) * invDir.z_;
    const float tzmax = (bounds[1-dirIsNeg[2]].z_ - ray.org_.z_) * invDir.z_;
    if ((tmin > tzmax) || (tzmin > tmax))
      return false;
    if (tzmin > tmin)
      tmin = tzmin;
    if (tzmax < tmax)
      tmax = tzmax;
    return (tmin < kINF) && (tmax > 0.0);
  }
 public:
  virtual ~BVHBase() {}
  virtual void create(const std::vector<RefTriangle*>& triangles) = 0;
  virtual bool intersect(const Ray &ray, Hitpoint* hitpoint) = 0;
};


#endif
