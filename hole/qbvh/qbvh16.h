#ifndef _QBVH_16_
#define _QBVH_16_

#include <immintrin.h>
#include "sse.h"
#include "qbvh_table.h"

#include "bvh_base.h"

class QBVH16 : public BVHBase {
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

	
  // 16個の三角形をパックする
  struct SIMDTrianglePack {
    __m512 x_[3];
    __m512 y_[3];
    __m512 z_[3];
    int idx_[16];
  };
	

  struct Children {
    union U{
      struct Node {
	unsigned flag_ : 1;
	unsigned index_: 31;
      } node_;
      struct Leaf {
	unsigned flag_ : 1;
	unsigned nPrimitives_: 5;
	unsigned index_: 26;
      } leaf_;

      unsigned int raw_;
    } u;
  };

  struct SIMDBVHNode{
    __m512 bboxes[2][3];//4 float min-max xyz
    Children children[4]; //4 children
    int axis_top;       //top axis
    int axis_left;      //left axis
    int axis_right;     //right axis
    int reserved;       //padding
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
  std::vector<SIMDBVHNode*> simd_nodes_;

  LinearBVHNode*nodes_;


  inline int test_AABB(
		       const __m512 bboxes[2][3],  //4boxes : min-max[2] of xyz[3] of boxes[4]
		       
		       const __m512 org[3],        //ray origin
		       const __m512 idir[3],       //ray inveresed direction
		       const int sign[3],          //ray xyz direction -> +:0,-:1
		       __m512 tmin, __m512 tmax    //ray range tmin-tmax 
		       )
  {
    // x coordinate
    tmin = _mm512_max_ps(
		      tmin,
		      _mm512_mul_ps(_mm512_sub_ps(bboxes[sign[0]][0],org[0]), idir[0])
		      );
    tmax = _mm512_min_ps(
		      tmax,
		      _mm512_mul_ps(_mm512_sub_ps(bboxes[1 - sign[0]][0], org[0]), idir[0])
		      );

    // y coordinate
    tmin = _mm512_max_ps(
		      tmin,
		      _mm512_mul_ps(_mm512_sub_ps(bboxes[sign[1]][1],org[1]), idir[1])
		      );
    tmax = _mm512_min_ps(
		      tmax,
		      _mm512_mul_ps(_mm512_sub_ps(bboxes[1 - sign[1]][1], org[1]), idir[1])
		      );

    // z coordinate
    tmin = _mm512_max_ps(
		      tmin,
		      _mm512_mul_ps(_mm512_sub_ps(bboxes[sign[2]][2],org[2]), idir[2])
		      );
    tmax = _mm512_min_ps(
		      tmax,
		      _mm512_mul_ps(_mm512_sub_ps(bboxes[1 - sign[2]][2], org[2]), idir[2])
		      );
    return _mm512_mask2int(_mm512_cmple_ps_mask(tmin, tmax));//tmin<tmaxとなれば交差
  }



  BVHBuildNode* build(std::vector<BVHPrimitiveInfo> &buildData, int start, int end, int *totalNodes) {
    (*totalNodes) ++;
    BVHBuildNode* node = new BVHBuildNode;
    BBox bbox;
    for (int i = start; i < end; ++i)
      bbox = Union(bbox, buildData[i].bounds_);

    int nPrimitives = end - start;
    if (nPrimitives <= 16) {
      // 葉
      int firstPrimOffset = ordered_triangles_.size();

      SIMDTrianglePack *simdt = (SIMDTrianglePack*)_aligned_malloc(sizeof(SIMDTrianglePack), 64);
			
      float x[16 * 3] __attribute__ ((aligned (64))) = {0};
      float y[16 * 3] __attribute__ ((aligned (64))) = {0};
      float z[16 * 3] __attribute__ ((aligned (64))) = {0};

      int cnt = 0;
      for (int i = start; i < end; ++i , ++cnt){
	const int idx = buildData[i].primitiveNumber_;
	ordered_triangles_.push_back(original_triangles_[idx]);
				
	int t = cnt % 16;

	simdt->idx_[t] = firstPrimOffset + cnt;
	x[t]     = original_triangles_[idx]->p_[0]->x_;
	x[16 + t] = original_triangles_[idx]->p_[1]->x_;
	x[32 + t] = original_triangles_[idx]->p_[2]->x_;

	y[t]     = original_triangles_[idx]->p_[0]->y_;
	y[16 + t] = original_triangles_[idx]->p_[1]->y_;
	y[32 + t] = original_triangles_[idx]->p_[2]->y_;

	z[t]     = original_triangles_[idx]->p_[0]->z_;
	z[16 + t] = original_triangles_[idx]->p_[1]->z_;
	z[32 + t] = original_triangles_[idx]->p_[2]->z_;
      }
      for (; cnt < 16; ++cnt) {
	simdt->idx_[cnt%16] = -1;
      }

      for (int i = 0; i < 3; ++i) {
	simdt->x_[i] = _mm512_load_ps(x + 16 * i);
	simdt->y_[i] = _mm512_load_ps(y + 16 * i);
	simdt->z_[i] = _mm512_load_ps(z + 16 * i);
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

  void collapse2QBVH(BVHBuildNode* node) {
    BVHBuildNode *lc = node->children_[0];
    BVHBuildNode *rc = node->children_[1];

    BVHBuildNode *c[4] = {0};
		
    SIMDBVHNode *n;
    n = (SIMDBVHNode*)_aligned_malloc(sizeof(SIMDBVHNode), 64);
    simd_nodes_.push_back(n);
    n->axis_top = node->splitAxis_;
    n->axis_left = n->axis_right = 0;

    if (lc != NULL) {
      n->axis_left = lc->splitAxis_;
      if (lc->nPrimitives_ == 0) {
	c[0] = lc->children_[0];
	c[1] = lc->children_[1];
      } else {
	c[0] = lc;
      }
    }
    if (rc != NULL) {
      n->axis_right = rc->splitAxis_;
      if (rc->nPrimitives_ == 0) {
	c[2] = rc->children_[0];
	c[3] = rc->children_[1];
      } else {
	c[2] = rc;
      }
    }
    float bboxes[2][3][16] __attribute__ ((aligned (64))) = {0};
    for (int j = 0; j < 3; ++j) {
      for (int k = 0; k < 4; ++k) {
	if (c[k] != NULL) {
	  /*
	    bboxes[0][j][k] = c[k]->bounds_.pmin_[j] - 0.2;
	    bboxes[1][j][k] = c[k]->bounds_.pmax_[j] + 0.2;
	  */
	  bboxes[0][j][k] = c[k]->bounds_.pmin_[j];
	  bboxes[1][j][k] = c[k]->bounds_.pmax_[j];
	}
      }
    }
    for(int m = 0; m < 2; m++){//minmax
      for(int j = 0; j < 3; j++){//xyz
	n->bboxes[m][j] = _mm512_load_ps(bboxes[m][j]);
      }
    }

    for (int i = 0; i < 4; ++i) {
      if (c[i] == NULL) {
	n->children[i].u.leaf_.flag_ = 1;
	n->children[i].u.leaf_.nPrimitives_ = 0;
	n->children[i].u.leaf_.index_ = 0;
      } else {
	if (c[i]->nPrimitives_ == 0) {
	  n->children[i].u.node_.flag_ = 0;
	  n->children[i].u.node_.index_= simd_nodes_.size();
	  collapse2QBVH(c[i]);
	} else {
	  // std::cout << c[i]->nPrimitives_ << " ";
	  n->children[i].u.leaf_.flag_ = 1;
	  n->children[i].u.leaf_.nPrimitives_ = c[i]->nPrimitives_;
	  n->children[i].u.leaf_.index_ = c[i]->simdTrisIdx_;
	}
      }
    }

    return;
  }
    
 public:
  void set_value(float *array, const float value, const int num) {
    for (int i = 0; i < num; ++i)
      array[i] = value;
  }

  QBVH16() {
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

    // collapse
    collapse2QBVH(build(build_data, 0, original_triangles_.size(), &totalNodes));
    /*
    BVHBuildNode* root = build(build_data, 0, original_triangles_.size(), &totalNodes);
		
    nodes_ = new LinearBVHNode[totalNodes];
    int offset = 0;
    flatten(root, &offset);
    */
  }

  bool intersect(const Ray &ray, Hitpoint* hitpoint) {
    __m512 zero_;
    __m512 one_;
    __m512 inf_;
    __m512 keps_;

    float one_f[16] __attribute__ ((aligned (64)));
    set_value(one_f, 1.0f, 16);
    float inf_f[16] __attribute__ ((aligned (64)));
    set_value(inf_f, kINF, 16);
    float zero_f[16] __attribute__ ((aligned (64)));
    set_value(zero_f, 0.0f, 16);
    float keps_f[16] __attribute__ ((aligned (64)));
    set_value(keps_f, kEPS, 16);
   
    zero_ = _mm512_load_ps(zero_f);
    one_ = _mm512_load_ps(one_f);
    inf_ = _mm512_load_ps(inf_f);
    keps_ = _mm512_load_ps(keps_f);
       
    const Vec invDir(1.0f / ray.dir_.x_, 1.0f / ray.dir_.y_, 1.0f / ray.dir_.z_);
    const int dirIsNeg[3] = {invDir.x_ < 0, invDir.y_ < 0, invDir.z_ < 0 };

    __m512 sseOrg[3];
    __m512 sseiDir[3];
    int sign[3];
    Vec idir(1.0f / ray.dir_.x_, 1.0f / ray.dir_.y_, 1.0f / ray.dir_.z_);
    float r_idir_x[16] __attribute__ ((aligned (64)));
    set_value(r_idir_x, idir.x_, 16);
    float r_idir_y[16] __attribute__ ((aligned (64)));
    set_value(r_idir_y, idir.y_, 16);
    float r_idir_z[16] __attribute__ ((aligned (64)));
    set_value(r_idir_z, idir.z_, 16);
		 			
    float r_org_x[16] __attribute__ ((aligned (64)));
    set_value(r_org_x, ray.org_.x_, 16);
    float r_org_y[16] __attribute__ ((aligned (64)));
    set_value(r_org_y, ray.org_.y_, 16);
    float r_org_z[16] __attribute__ ((aligned (64)));
    set_value(r_org_z, ray.org_.z_, 16);

    float r_dir_x[16] __attribute__ ((aligned (64)));
    set_value(r_dir_x, ray.dir_.x_, 16);
    float r_dir_y[16] __attribute__ ((aligned (64)));
    set_value(r_dir_y, ray.dir_.y_, 16);
    float r_dir_z[16] __attribute__ ((aligned (64)));
    set_value(r_dir_z, ray.dir_.z_, 16);

    __m512 dir_x = _mm512_load_ps(r_dir_x);
    __m512 dir_y = _mm512_load_ps(r_dir_y);
    __m512 dir_z = _mm512_load_ps(r_dir_z);

    sseOrg[0] = _mm512_load_ps(r_org_x);
    sseOrg[1] = _mm512_load_ps(r_org_y);
    sseOrg[2] = _mm512_load_ps(r_org_z);

    sseiDir[0] = _mm512_load_ps(r_idir_x);
    sseiDir[1] = _mm512_load_ps(r_idir_y);
    sseiDir[2] = _mm512_load_ps(r_idir_z);
		
    sign[0] = idir[0] < 0;
    sign[1] = idir[1] < 0;
    sign[2] = idir[2] < 0;

    const SIMDBVHNode* nodes = simd_nodes_[0];
    
    Children nodeStack[40];
    int todoNode = 0;
    nodeStack[0].u.raw_ = 0;

    bool hit = false;
    int triangle_index = -1;

    while (true) {
      Children item = nodeStack[todoNode];
      todoNode--;//pop stack

      if(item.u.node_.flag_ == 0){
	const SIMDBVHNode& node = *(simd_nodes_[item.u.node_.index_]);
	float now_distance_f[16] __attribute__ ((aligned (64)));
	set_value(now_distance_f, hitpoint->distance_, 16);
	__m512 now_distance = _mm512_load_ps(now_distance_f);
	const int HitMask = test_AABB(node.bboxes, sseOrg, sseiDir, sign, zero_, now_distance);

	if (HitMask) {
	  const int nodeIdx = (sign[node.axis_top] << 2) | (sign[node.axis_left] << 1) | (sign[node.axis_right]);
	  int bboxOrder = OrderTable[HitMask * 8 + nodeIdx];
					
	  for (int i = 0; i < 4; ++i) {
	    if (bboxOrder & 0x4)
	      break;
	    ++todoNode;
	    nodeStack[todoNode] = node.children[bboxOrder & 0x3];

	    bboxOrder >>= 4;
	  }
	}


      } else {
	int no_hit_f __attribute__ ((aligned (64)));
	float t_f[16] __attribute__ ((aligned (64)));
	int nohitmask;
	SIMDTrianglePack *s = simd_triangles_[item.u.leaf_.index_];
		
	__m512 e1_x = _mm512_sub_ps(s->x_[1], s->x_[0]);
	__m512 e1_y = _mm512_sub_ps(s->y_[1], s->y_[0]);
	__m512 e1_z = _mm512_sub_ps(s->z_[1], s->z_[0]);
					
	__m512 e2_x = _mm512_sub_ps(s->x_[2], s->x_[0]);
	__m512 e2_y = _mm512_sub_ps(s->y_[2], s->y_[0]);
	__m512 e2_z = _mm512_sub_ps(s->z_[2], s->z_[0]);

	__m512 s1_x = _mm512_sub_ps(_mm512_mul_ps(dir_y, e2_z), _mm512_mul_ps(dir_z, e2_y));
	__m512 s1_y = _mm512_sub_ps(_mm512_mul_ps(dir_z, e2_x), _mm512_mul_ps(dir_x, e2_z));
	__m512 s1_z = _mm512_sub_ps(_mm512_mul_ps(dir_x, e2_y), _mm512_mul_ps(dir_y, e2_x));

	__m512 divisor = _mm512_add_ps(_mm512_add_ps(_mm512_mul_ps(s1_x, e1_x), _mm512_mul_ps(s1_y, e1_y)), _mm512_mul_ps(s1_z, e1_z));
	__mmask16 no_hit  = _mm512_cmpeq_ps_mask(divisor, zero_);	

	__m512 invDivisor = _mm512_rcp23_ps(divisor);
	//	  __m512 invDivisor = _mm512_div_ps(one_, divisor);

	__m512 d_x = _mm512_sub_ps(sseOrg[0], s->x_[0]);
	__m512 d_y = _mm512_sub_ps(sseOrg[1], s->y_[0]);
	__m512 d_z = _mm512_sub_ps(sseOrg[2], s->z_[0]);

	__m512 b1 = _mm512_mul_ps(_mm512_add_ps(_mm512_add_ps(_mm512_mul_ps(d_x, s1_x), _mm512_mul_ps(d_y, s1_y)), _mm512_mul_ps(d_z, s1_z)),
				  invDivisor);
	no_hit = _mm512_kor(no_hit, _mm512_kor(_mm512_cmplt_ps_mask(b1, zero_), _mm512_cmplt_ps_mask(one_, b1)));
					
	__m512 s2_x = _mm512_sub_ps(_mm512_mul_ps(d_y, e1_z), _mm512_mul_ps(d_z, e1_y));
	__m512 s2_y = _mm512_sub_ps(_mm512_mul_ps(d_z, e1_x), _mm512_mul_ps(d_x, e1_z));
	__m512 s2_z = _mm512_sub_ps(_mm512_mul_ps(d_x, e1_y), _mm512_mul_ps(d_y, e1_x));

	__m512 b2 = _mm512_mul_ps(_mm512_add_ps(_mm512_add_ps(_mm512_mul_ps(dir_x, s2_x), _mm512_mul_ps(dir_y, s2_y)), _mm512_mul_ps(dir_z, s2_z)),
				  invDivisor);
	no_hit = _mm512_kor(no_hit, _mm512_kor(_mm512_cmplt_ps_mask(b2, zero_), _mm512_cmplt_ps_mask(one_, _mm512_add_ps(b1, b2))));

	__m512 t = _mm512_mul_ps(_mm512_add_ps(_mm512_add_ps(_mm512_mul_ps(e2_x, s2_x), _mm512_mul_ps(e2_y, s2_y)), _mm512_mul_ps(e2_z, s2_z)),
				 invDivisor);
				
	no_hit = _mm512_kor(no_hit, _mm512_cmplt_ps_mask(t, keps_));
				
	nohitmask = _mm512_mask2int(no_hit);
	_mm512_store_ps(t_f, t);
				
	for (int i = 0; i < 16; ++i) {
	  if ((nohitmask & (1 << i)) == 0 && hitpoint->distance_ > t_f[i]) {
	    hit = true;
	    triangle_index = s->idx_[i];
	    hitpoint->distance_ = t_f[i];
	  }
	}
      }

      if (todoNode < 0)
	break;
    }
    
    if (hit) {
      if (triangle_index >= ordered_triangles_.size()) {
	std::cout << "HOGE";
	return false;
      }
      RefTriangle *t = ordered_triangles_[triangle_index];
      hitpoint->position_ = ray.org_+ hitpoint->distance_ * ray.dir_;
      hitpoint->normal_   = normalize(cross((*t->p_[2]) - (*t->p_[0]), (*t->p_[1]) - (*t->p_[0])));
    }

    return hit;
  }

  virtual ~QBVH16() {
  }
};

#endif
