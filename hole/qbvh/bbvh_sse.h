#ifndef _BBVH_SSE_
#define _BBVH_SSE_

#include <smmintrin.h>
#include "sse.h"

#include "bvh_base.h"

class BBVHSSE : public BVHBase {
 private:
  struct BVHBuildNode  {
    BBox bounds_;
    BVHBuildNode* children_[2];
    int splitAxis_, firstPrimOffset_, nPrimitives_;
    int simdTrisIdx_;

    BVHBuildNode() {
    }
		
    void init_leaf(int first, int n, const BBox& b, const int simdTrisIdx) {
      firstPrimOffset_ = first;
      nPrimitives_ = n;
      bounds_ = b;
      simdTrisIdx_ = simdTrisIdx;
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
    int simdTrisIdx_;
    // char pad[2];
  };

	
  // 四つの三角形をパックする
  struct SIMDTrianglePack {
    __m128 x_[3];
    __m128 y_[3];
    __m128 z_[3];
    int idx_[4];
  };
	
  struct SIMDLinearBVHNode {
    BBox bounds_;
    union {
      int primitiveOffset_;
      int secondChildOffset_;
    };
    char nPrimitives_;
    char axis_;
    int simdTrisIdx_;
  };

  std::vector<RefTriangle*> ordered_triangles_;
  std::vector<RefTriangle*> original_triangles_;
  std::vector<SIMDTrianglePack*> simd_triangles_;
  LinearBVHNode*nodes_;

  BVHBuildNode* build(std::vector<BVHPrimitiveInfo> &buildData, int start, int end, int *totalNodes) {
    (*totalNodes) ++;
    BVHBuildNode* node = new BVHBuildNode;
    BBox bbox;
    for (int i = start; i < end; ++i)
      bbox = Union(bbox, buildData[i].bounds_);

    int nPrimitives = end - start;
    if (nPrimitives <= 4) {
      // 葉
      int firstPrimOffset = ordered_triangles_.size();

      SIMDTrianglePack *simdt = (SIMDTrianglePack*)_aligned_malloc(sizeof(SIMDTrianglePack), 16);
			
      float x[4 * 3] __attribute__ ((aligned (16))) = {0};
      float y[4 * 3] __attribute__ ((aligned (16))) = {0};
      float z[4 * 3] __attribute__ ((aligned (16))) = {0};

      int cnt = 0;
      for (int i = start; i < end; ++i , ++cnt){
	const int idx = buildData[i].primitiveNumber_;
	ordered_triangles_.push_back(original_triangles_[idx]);
				
	int t = cnt % 4;

	simdt->idx_[t] = firstPrimOffset + cnt;
	x[t]     = original_triangles_[idx]->p_[0]->x_;
	x[4 + t] = original_triangles_[idx]->p_[1]->x_;
	x[8 + t] = original_triangles_[idx]->p_[2]->x_;

	y[t]     = original_triangles_[idx]->p_[0]->y_;
	y[4 + t] = original_triangles_[idx]->p_[1]->y_;
	y[8 + t] = original_triangles_[idx]->p_[2]->y_;

	z[t]     = original_triangles_[idx]->p_[0]->z_;
	z[4 + t] = original_triangles_[idx]->p_[1]->z_;
	z[8 + t] = original_triangles_[idx]->p_[2]->z_;
      }
      for (; cnt < 4; ++cnt) {
	simdt->idx_[cnt%4] = -1;
      }

      for (int i = 0; i < 3; ++i) {
	simdt->x_[i] = _mm_load_ps(x + 4 * i);
	simdt->y_[i] = _mm_load_ps(y + 4 * i);
	simdt->z_[i] = _mm_load_ps(z + 4 * i);
      }

      simd_triangles_.push_back(simdt);
      node->init_leaf(firstPrimOffset, nPrimitives, bbox, simd_triangles_.size() - 1);
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
      linearNode->simdTrisIdx_ = node->simdTrisIdx_;
    } else {
      linearNode->axis_ = node->splitAxis_;
      linearNode->nPrimitives_ = 0;
      flatten(node->children_[0], offset);
      linearNode->secondChildOffset_ = flatten(node->children_[1], offset);
    }

    return myOffset;
  }
	
  __m128 zero_;
  __m128 one_;
  __m128 inf_;
  __m128 keps_;
 public:
  BBVHSSE() {
    const float one_f[4] __attribute__ ((aligned (16))) = {1.0f, 1.0f, 1.0f, 1.0f};
    const float inf_f[4] __attribute__ ((aligned (16))) = {kINF, kINF, kINF, kINF};
    const float zero_f[4] __attribute__ ((aligned (16)))= {0.0f, 0.0f, 0.0f, 0.0f};
    const float keps_f[4] __attribute__ ((aligned (16)))= {kEPS, kEPS, kEPS, kEPS};
		
    zero_ = _mm_load_ps(zero_f);
    one_ = _mm_load_ps(one_f);
    inf_ = _mm_load_ps(inf_f);
    keps_ = _mm_load_ps(keps_f);
  }

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

    __m128 sseOrg[3];
    __m128 sseiDir[3];
    int sign[3];
    Vec idir(1.0f / ray.dir_.x_, 1.0f / ray.dir_.y_, 1.0f / ray.dir_.z_);
    float r_idir_x[4] __attribute__ ((aligned (16))) = {idir.x_, idir.x_, idir.x_, idir.x_};
    float r_idir_y[4] __attribute__ ((aligned (16))) = {idir.y_, idir.y_, idir.y_, idir.y_};
    float r_idir_z[4] __attribute__ ((aligned (16))) = {idir.z_, idir.z_, idir.z_, idir.z_};
		 			
    float r_org_x[4] __attribute__ ((aligned (16)))= {ray.org_.x_, ray.org_.x_, ray.org_.x_, ray.org_.x_};
    float r_org_y[4] __attribute__ ((aligned (16)))= {ray.org_.y_, ray.org_.y_, ray.org_.y_, ray.org_.y_};
    float r_org_z[4] __attribute__ ((aligned (16)))= {ray.org_.z_, ray.org_.z_, ray.org_.z_, ray.org_.z_};

    float r_dir_x[4] __attribute__ ((aligned (16)))= {ray.dir_.x_, ray.dir_.x_, ray.dir_.x_, ray.dir_.x_};
    float r_dir_y[4] __attribute__ ((aligned (16)))= {ray.dir_.y_, ray.dir_.y_, ray.dir_.y_, ray.dir_.y_};
    float r_dir_z[4] __attribute__ ((aligned (16)))= {ray.dir_.z_, ray.dir_.z_, ray.dir_.z_, ray.dir_.z_};
					
    __m128 dir_x = _mm_load_ps(r_dir_x);
    __m128 dir_y = _mm_load_ps(r_dir_y);
    __m128 dir_z = _mm_load_ps(r_dir_z);

    sseOrg[0] = _mm_load_ps(r_org_x);
    sseOrg[1] = _mm_load_ps(r_org_y);
    sseOrg[2] = _mm_load_ps(r_org_z);

    sseiDir[0] = _mm_load_ps(r_idir_x);
    sseiDir[1] = _mm_load_ps(r_idir_y);
    sseiDir[2] = _mm_load_ps(r_idir_z);
		
    sign[0] = idir[0] < 0;
    sign[1] = idir[1] < 0;
    sign[2] = idir[2] < 0;

    int todoOffset = 0, nodeNum = 0;
    int todo[256];

    bool hit = false;
    int triangle_index = -1;

    while (true) {
      const LinearBVHNode *node = &nodes_[nodeNum];

      if (intersect_bbox(node->bounds_, ray, invDir, dirIsNeg)) {

	if (node->nPrimitives_ > 0) {
	  /*
	    for (int j = 0; j < node->nPrimitives_; ++j) {
	    Hitpoint tmp_hitpoint;
	    if (ordered_triangles_[node->primitiveOffset_ + j]->intersect(ray, &tmp_hitpoint)) {
	    if (hitpoint->distance_ > tmp_hitpoint.distance_) {
	    hit = true;
	    *hitpoint = tmp_hitpoint;
	    }
	    }
	    }*/

	  // SSE!
					
	  float no_hit_f[4] __attribute__ ((aligned (16)));
	  float t_f[4] __attribute__ ((aligned (16)));
	  int nohitmask;
	  SIMDTrianglePack *s = simd_triangles_[node->simdTrisIdx_];

				
	  __m128 e1_x = _mm_sub_ps(s->x_[1], s->x_[0]);
	  __m128 e1_y = _mm_sub_ps(s->y_[1], s->y_[0]);
	  __m128 e1_z = _mm_sub_ps(s->z_[1], s->z_[0]);
					
	  __m128 e2_x = _mm_sub_ps(s->x_[2], s->x_[0]);
	  __m128 e2_y = _mm_sub_ps(s->y_[2], s->y_[0]);
	  __m128 e2_z = _mm_sub_ps(s->z_[2], s->z_[0]);

	  __m128 s1_x = _mm_sub_ps(_mm_mul_ps(dir_y, e2_z), _mm_mul_ps(dir_z, e2_y));
	  __m128 s1_y = _mm_sub_ps(_mm_mul_ps(dir_z, e2_x), _mm_mul_ps(dir_x, e2_z));
	  __m128 s1_z = _mm_sub_ps(_mm_mul_ps(dir_x, e2_y), _mm_mul_ps(dir_y, e2_x));

	  __m128 divisor = _mm_add_ps(_mm_add_ps(_mm_mul_ps(s1_x, e1_x), _mm_mul_ps(s1_y, e1_y)), _mm_mul_ps(s1_z, e1_z));
	  __m128 no_hit  = _mm_cmpeq_ps(divisor, zero_);	

	  __m128 invDivisor = _mm_rcp_ps(divisor);

	  __m128 d_x = _mm_sub_ps(sseOrg[0], s->x_[0]);
	  __m128 d_y = _mm_sub_ps(sseOrg[1], s->y_[0]);
	  __m128 d_z = _mm_sub_ps(sseOrg[2], s->z_[0]);

	  __m128 b1 = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(d_x, s1_x), _mm_mul_ps(d_y, s1_y)), _mm_mul_ps(d_z, s1_z)),
				 invDivisor);
	  no_hit = _mm_or_ps(no_hit, _mm_or_ps(_mm_cmplt_ps(b1, zero_), _mm_cmpgt_ps(b1, one_)));
					
	  __m128 s2_x = _mm_sub_ps(_mm_mul_ps(d_y, e1_z), _mm_mul_ps(d_z, e1_y));
	  __m128 s2_y = _mm_sub_ps(_mm_mul_ps(d_z, e1_x), _mm_mul_ps(d_x, e1_z));
	  __m128 s2_z = _mm_sub_ps(_mm_mul_ps(d_x, e1_y), _mm_mul_ps(d_y, e1_x));

	  __m128 b2 = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(dir_x, s2_x), _mm_mul_ps(dir_y, s2_y)), _mm_mul_ps(dir_z, s2_z)),
				 invDivisor);
	  no_hit = _mm_or_ps(no_hit, _mm_or_ps(_mm_cmplt_ps(b2, zero_), _mm_cmpgt_ps(_mm_add_ps(b1, b2), one_)));

	  __m128 t = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(e2_x, s2_x), _mm_mul_ps(e2_y, s2_y)), _mm_mul_ps(e2_z, s2_z)),
				invDivisor);
				
	  no_hit = _mm_or_ps(no_hit, _mm_cmplt_ps(t, keps_));
				
	  nohitmask = _mm_movemask_ps(no_hit);
	  _mm_store_ps(t_f, t);
				
	  for (int i = 0; i < 4; ++i) {
	    if ((nohitmask & (1 << i)) == 0 && hitpoint->distance_ > t_f[i]) {
	      hit = true;
	      triangle_index = s->idx_[i];
	      hitpoint->distance_ = t_f[i];
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
		
    if (hit) {
      RefTriangle *t = ordered_triangles_[triangle_index];
      hitpoint->position_ = ray.org_+ hitpoint->distance_ * ray.dir_;
      hitpoint->normal_   = normalize(cross((*t->p_[2]) - (*t->p_[0]), (*t->p_[1]) - (*t->p_[0])));
    }

    return hit;
  }

  virtual ~BBVHSSE() {
  }
};

#endif
