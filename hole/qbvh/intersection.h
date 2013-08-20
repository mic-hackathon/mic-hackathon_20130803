#ifndef _INTERSECTION_H_
#define _INTERSECTION_H_

#include "vec.h"
#include "constant.h"

struct Hitpoint {
	float distance_;
	Vec normal_;
	Vec position_;

	
	int triangle_id;
	int object_id;

	Hitpoint() : distance_(kINF), normal_(), position_(), triangle_id(-1), object_id(-1) {}
};


#endif
