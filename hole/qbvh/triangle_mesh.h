#ifndef _TRIANGLE_MESH_
#define _TRIANGLE_MESH_

#include <iostream>
#include <vector>

#include "triangle.h"

#include "bbvh.h"

#ifdef __MIC__
#include "bbvh16.h"
#include "nbvh16.h"
#include "qbvh16.h"
#endif
//#include "bbvh_sse.h"
#include "qbvh_sse.h"

enum BVHType {
  BVHType_BBVH,
  BVHType_BBVHSSE,
  BVHType_QBVHSSE,

  BVHType_BBVH16,
  BVHType_QBVH16,
  BVHType_NBVH16,
};

class TriangleMesh {	
 private:
  std::vector<Vec> vertex;
  std::vector<Vec> face;

  BBox bbox;

  // BBVH bvh;
  //BBVHSSE bvh;
  //  QBVHSSE bvh;
  BVHBase *bvh;
 public:
  std::vector<RefTriangle*> triangles;
  TriangleMesh(std::string filename, double scale, Vec origin, BVHType type) {
    FILE* fp = fopen(filename.c_str(), "rt");
    if (fp == NULL) {
      std::cout << "error: " << filename << std::endl;
      return;
    }

    if (fp != NULL) {
      std::cout << "Load obj file: " << filename <<  std::endl;
      char buf[1024];
      while (fgets(buf, 1024, fp) != NULL) {
	char prefix = buf[0];
	switch (prefix) {
	case '#':
	  continue;
	case 'o':
	  continue;
	case 'v':
	  {
	    float v1, v2, v3;
	    sscanf(buf, "v %f %f %f", &v1, &v2, &v3);

	    Vec pt = Vec(v1, v2, v3);
	    const Vec v = origin + scale * pt;
	    vertex.push_back(v);
	  }
	  break;
	case 'f':
	  {
	    int f1, f2, f3;
	    sscanf(buf, "f %d %d %d", &f1, &f2, &f3);
	    face.push_back(Vec(f3-1, f2-1, f1-1));
	  }
	  break;

	default:
	  break;
	}
      }

      fclose(fp);
    }


    std::cout << "vertex: " << vertex.size() << std::endl;
    std::cout << "face: " << face.size() << std::endl;

    for (int i = 0; i < face.size(); i ++) {
      triangles.push_back(new RefTriangle(&vertex[face[i][0]], &vertex[face[i][1]], &vertex[face[i][2]]));
    }

    for (int i = 0; i < vertex.size(); i ++) {
      bbox = Union(bbox, vertex[i]);
    }

    std::cout << bbox.pmin_.x_ << " " << bbox.pmin_.y_ << " " << bbox.pmin_.z_ << std::endl;
    std::cout << bbox.pmax_.x_ << " " << bbox.pmax_.y_ << " " << bbox.pmax_.z_ << std::endl;

    std::cout << "======================== ";
    switch (type) {
    case BVHType_BBVH:
      std::cout << "BBVH" << std::endl;
      bvh = new BBVH();
      break;
    case BVHType_BBVHSSE:
      std::cout << "BBVH SSE" << std::endl;
      //bvh = new BBVHSSE();
      break;
    case BVHType_QBVHSSE:
      std::cout << "QBVH SSE" << std::endl;
      bvh = new QBVHSSE();
      break;

    case BVHType_BBVH16:
      #ifdef __MIC__
      std::cout << "BBVH16" << std::endl;
      bvh = new BBVH16();
      #endif
      break;

    case BVHType_QBVH16:
      #ifdef __MIC__
      std::cout << "QBVH16" << std::endl;
      bvh = new QBVH16();
      #endif
      break;

    case BVHType_NBVH16:
      #ifdef __MIC__
      std::cout << "NBVH16" << std::endl;
      bvh = new NBVH16();
      #endif
      break;
    }
		
    std::cout << "create BVH" << std::endl;
    bvh->create(triangles);
    std::cout << "done BVH" << std::endl << std::endl;
  }

  virtual ~TriangleMesh() {
    for (int i = 0; i < triangles.size(); i ++)
      delete triangles[i];
  }

  BBox ObjectBound() const {
    return bbox;
  }

  bool intersect(const Ray &ray, Hitpoint* hitpoint) {
    bool check = bbox.check_intersect(ray, NULL, NULL);
    if (check) {
      const bool ret = bvh->intersect(ray, hitpoint);
      return ret;
    }

    return false;
  }
	
};

#endif
