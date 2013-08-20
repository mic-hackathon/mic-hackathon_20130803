#ifndef _BBVH_
#define _BBVH_

#include "bvh_base.h"

class BBVH : public BVHBase {
 private:
  struct BVHBuildNode  {
    BBox bounds_;
    BVHBuildNode* children_[2];
    int splitAxis_, firstPrimOffset_, nPrimitives_;

    BVHBuildNode() {
    }

    void init_leaf(int first, int n, const BBox& b) {
      firstPrimOffset_ = first;
      nPrimitives_ = n;
      bounds_ = b;
      splitAxis_ = 0;
    }

    void init_interior(int axis, BVHBuildNode* c0, BVHBuildNode* c1) {
      children_[0] = c0;
      children_[1] = c1;
      bounds_ = Union(c0->bounds_, c1->bounds_);
      splitAxis_ = axis;
      firstPrimOffset_ = -1;
      nPrimitives_ = 0;
    }
  };
	
  struct LinearBVHNode {
    BBox bounds_;
    union {
      int primitiveOffset_;
      int secondChildOffset_;
    };
    char nPrimitives_;
    char axis_;
    // char pad[2];
  };


  std::vector<RefTriangle*> ordered_triangles_;
  std::vector<RefTriangle*> original_triangles_;
  LinearBVHNode*nodes_;

  BVHBuildNode* build(std::vector<BVHPrimitiveInfo> &buildData, int start, int end, int *totalNodes) {
    (*totalNodes) ++;
    BVHBuildNode* node = new BVHBuildNode;
    BBox bbox;
    for (int i = start; i < end; ++i)
      bbox = Union(bbox, buildData[i].bounds_);

    int nPrimitives = end - start;
    if (nPrimitives <= 4) {
      // —t
      int firstPrimOffset = ordered_triangles_.size();
      for (int i = start; i < end; ++i ){
	int primNum = buildData[i].primitiveNumber_;
	ordered_triangles_.push_back(original_triangles_[primNum]);
      }
      node->init_leaf(firstPrimOffset, nPrimitives, bbox);
    } else {
      int mid;
      int dim;
			
      split(buildData, bbox, start, end, &mid, &dim);

      node->init_interior(dim,
			  build(buildData, start, mid, totalNodes),
			  build(buildData, mid, end, totalNodes));
    }

    return node;
  }
	
  int flatten(BVHBuildNode* node, int *offset) {
    LinearBVHNode* linearNode = &nodes_[*offset];
    linearNode->bounds_ = node->bounds_;
    int myOffset = (*offset) ++;
    if (node->nPrimitives_ > 0) {
      linearNode->primitiveOffset_ = node->firstPrimOffset_;
      linearNode->nPrimitives_ = node->nPrimitives_;
    } else {
      linearNode->axis_ = node->splitAxis_;
      linearNode->nPrimitives_ = 0;
      flatten(node->children_[0], offset);
      linearNode->secondChildOffset_ = flatten(node->children_[1], offset);
    }

    return myOffset;
  }

 public:
  BBVH() {}

  void create(const std::vector<RefTriangle*>& triangles) {
    original_triangles_ = triangles;
    ordered_triangles_.clear();

    std::vector<BVHPrimitiveInfo> build_data;
    for (int i = 0; i < original_triangles_.size(); i ++) {
      BBox b = original_triangles_[i]->ObjectBound();
      build_data.push_back(BVHPrimitiveInfo(i, b));
    }
    int totalNodes = 0;
    ordered_triangles_.reserve(original_triangles_.size());

    BVHBuildNode* root = build(build_data, 0, original_triangles_.size(), &totalNodes);
		
    nodes_ = new LinearBVHNode[totalNodes];
    int offset = 0;
    flatten(root, &offset);
  }

  bool intersect(const Ray &ray, Hitpoint* hitpoint) {
    const Vec invDir(1.0f / ray.dir_.x_, 1.0f / ray.dir_.y_, 1.0f / ray.dir_.z_);
    const int dirIsNeg[3] = {invDir.x_ < 0, invDir.y_ < 0, invDir.z_ < 0 };

    int todoOffset = 0, nodeNum = 0;
    int todo[256];

    bool hit = false;

    while (true) {
      const LinearBVHNode *node = &nodes_[nodeNum];

      if (intersect_bbox(node->bounds_, ray, invDir, dirIsNeg)) {

	if (node->nPrimitives_ > 0) {
	  for (int j = 0; j < node->nPrimitives_; ++j) {
	    Hitpoint tmp_hitpoint;
	    if (ordered_triangles_[node->primitiveOffset_ + j]->intersect(ray, &tmp_hitpoint)) {
	      if (hitpoint->distance_ > tmp_hitpoint.distance_) {
		hit = true;
		*hitpoint = tmp_hitpoint;
	      }
	    }
	  }

	  if (todoOffset == 0) break;
	  nodeNum = todo[--todoOffset];
	} else {
	  if (dirIsNeg[node->axis_]) {
	    todo[todoOffset ++] = nodeNum + 1;
	    nodeNum = node->secondChildOffset_;
	  } else {
	    todo[todoOffset ++] = node->secondChildOffset_;
	    nodeNum = nodeNum + 1;
	  }
	}
      } else {
	if (todoOffset == 0) break;
	nodeNum = todo[--todoOffset];
      }
    }
    return hit;
  }

  virtual ~BBVH() {
  }
};

#endif
