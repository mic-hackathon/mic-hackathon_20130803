#include <iostream>
#include <omp.h>
#include <sys/time.h>

#include "vec.h"
#include "hdr.h"
#include "triangle_mesh.h"
#include "random.h"
#include "ibl.h"

// 時間計測関数
double time_diff(timeval *begin, timeval *end) {
  return (double)(end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec) / 1000000.0;
}
class Timer {
private:
  timeval t0_, t1_;
  std::string name_;
public:
  Timer(std::string name) : name_(name) {
    gettimeofday(&t0_, NULL);
  }

  ~Timer() {
    gettimeofday(&t1_, NULL);
    std::cout << name_ << ": " << time_diff(&t0_, &t1_) << " sec" << std::endl;
  }
};

IBL ibl("grace_probe_high.hdr", true, 64);

TriangleMesh* mesh;
struct Scene {
  TriangleMesh *mesh_;
  Vec color_;
};

std::vector<Scene> scene;

bool intersect_scene(const Ray &ray, Hitpoint *hitpoint, int *id) {
  *id = -1;

  //  std::cerr << "_";
  // シーンと交差判定
  for (int i = 0; i < scene.size(); ++i) {
    if (scene[i].mesh_->intersect(ray, hitpoint))
      *id = i;
  }

  if (*id == -1)
    return false;

  return true;
}

Color radiance(const Ray &ray, Random *rnd, const int depth) {
  if (depth == 5)
    return Color();

  Hitpoint hitpoint;
  int lastid;
 
  if (!intersect_scene(ray, &hitpoint, &lastid)) {
    return Color(1, 1, 1);

    /*
      const float theta = acos(ray.dir_.y_);
      float phi = acos(ray.dir_.x_ / sqrt(ray.dir_.x_ * ray.dir_.x_ + ray.dir_.z_ * ray.dir_.z_));
      if (ray.dir_.z_ < 0.0)
      phi = 2.0 * kPI - phi;
		
      return ibl.sample(phi, theta);
    */
  }


  const Vec orienting_normal = dot(hitpoint.normal_ , ray.dir_) < 0.0 ? hitpoint.normal_: (-1.0 * hitpoint.normal_); // 交差位置の法線（物体からのレイの入出を考慮）

  float russian_roulette_probability = 1.0f;

  Color incoming_radiance;
  Color weight = 1.0;

  // orienting_normalの方向を基準とした正規直交基底(w, u, v)を作る。この基底に対する半球内で次のレイを飛ばす。
  Vec w, u, v;
  w = orienting_normal;
  if (fabs(w.x_) > kEPS) // ベクトルwと直交するベクトルを作る。w.xが0に近い場合とそうでない場合とで使うベクトルを変える。
    u = normalize(cross(Vec(0.0, 1.0, 0.0), w));
  else
    u = normalize(cross(Vec(1.0, 0.0, 0.0), w));
  v = cross(w, u);
  // コサイン項を使った重点的サンプリング
  const float r1 = 2 * kPI * rnd->next01();
  const float r2 = rnd->next01(), r2s = sqrt(r2);
  Vec dir = normalize((
		       u * cos(r1) * r2s +
		       v * sin(r1) * r2s +
		       w * sqrt(1.0f - r2)));

  incoming_radiance = radiance(Ray(hitpoint.position_, dir), rnd, depth+1);
  weight = scene[lastid].color_ / russian_roulette_probability;
  return multiply(weight, incoming_radiance);
}

static Vec sun_dir = normalize(Vec(-0.5, 0.4, 0.7));
static Color sun_color(1, 1, 1);

Color direct_radiance(const Ray &ray, Random *rnd, const int depth) {
  if (depth == 5)
    return Color();

  Hitpoint hitpoint;
  int lastid = -1;

  if (!intersect_scene(ray, &hitpoint, &lastid)) {
    return Color(1, 1, 1);
  }
  const Vec orienting_normal = dot(hitpoint.normal_ , ray.dir_) < 0.0 ? hitpoint.normal_: (-1.0 * hitpoint.normal_); // 交差位置の法線（物体からのレイの入出を考慮）

  Color weight = 1.0;

  Hitpoint shadow_hitpoint;
  int shadow_id = -1;
  if (intersect_scene(Ray(hitpoint.position_, sun_dir), &shadow_hitpoint, &shadow_id)) {
    return Color(0, 0, 0);
  }

  weight = scene[lastid].color_;
  return multiply(weight, std::max(0.0f, dot(sun_dir, orienting_normal)) * sun_color);
}



int main(int argc, char **argv) {
  if (argc <= 9) {
    std::cout << "threadnum width height supersamples obj_filename r g b bvh_type{bbvh|bbvh_sse|qbvh_sse}" << std::endl;
    return 0;
  }
  const int num_thread = atoi(argv[1]);

  const int width = atoi(argv[2]); 
  const int height = atoi(argv[3]);
  const int supersamples = atoi(argv[4]);

  const std::string input_obj = argv[5];
  const Color input_color(atof(argv[6]), atof(argv[7]), atof(argv[8]));

  BVHType bvh_type = BVHType_BBVH;

  if (strcmp(argv[9], "bbvh") == 0)
    bvh_type = BVHType_BBVH;
  else if (strcmp(argv[9], "bbvh_sse") == 0)
    bvh_type = BVHType_BBVHSSE;
  else if (strcmp(argv[9], "qbvh_sse") == 0)
    bvh_type = BVHType_QBVHSSE;
  else if (strcmp(argv[9], "bbvh16") == 0)
    bvh_type = BVHType_BBVH16;
  else if (strcmp(argv[9], "qbvh16") == 0)
    bvh_type = BVHType_QBVH16;
  else if (strcmp(argv[9], "nbvh16") == 0)
    bvh_type = BVHType_NBVH16;

  HDRImage hdr(width, height);

  // カメラ位置
  const float angle = 0.05f * 2.0f * kPI;
  const Vec camera_position = Vec(50.0f + 0.9 * sin(angle), 0.55f, 50.0f + 0.9 * cos(angle));
  const Vec camera_dir = normalize(Vec(50.0f, 0.25f, 50.0f) - camera_position);
  const Vec camera_up = Vec(0.0f, 1.0f, 0.0f);
	
  // ワールド座標系でのスクリーンの大きさ
  const float screen_width = 30.0f * width / height;
  const float screen_height= 30.0f;
  // スクリーンまでの距離
  const float screen_dist  = 30.0f;
  // スクリーンを張るベクトル
  const Vec screen_x = normalize(cross(camera_dir, camera_up)) * screen_width;
  const Vec screen_y = normalize(cross(screen_x, camera_dir)) * screen_height;
  const Vec screen_center = camera_position + camera_dir * screen_dist;
	

  {
    Timer timer(std::string("===== load"));
    Scene sc;
    
    sc.mesh_ = new TriangleMesh(input_obj, 0.16f, Vec(50.0f, 0.0f, 50.0f), bvh_type);
    sc.color_ = input_color;
    scene.push_back(sc);

    sc.mesh_ = new TriangleMesh("cube.obj", 0.025f, Vec(50.0f, 0.0f, 50.0f), bvh_type);
    sc.color_ = Color(0.1f, 0.1f, 0.1f);
    scene.push_back(sc);
  }
  {
    Timer timer(std::string("===== rendering"));

    // 普通のパストレ＋環境マップ
    int fg = 0;
    
    // OpenMP

    omp_set_num_threads(num_thread);
#pragma omp parallel for
    for (int y = 0; y < height; y ++) {
      //		std::cerr << "Rendering (y = " << y << ") " << (100.0 * y / (height - 1)) << "%" << std::endl;

      if (y % 10 == 0)
	std::cerr << y << " ";

      Random rnd(y + 1);
      for (int x = 0; x < width; x ++) {
	// supersamples x supersamples のスーパーサンプリング
	for (int sy = 0; sy < supersamples; sy ++) {
	  for (int sx = 0; sx < supersamples; sx ++) {

	    Color accumulated_radiance = Color();
	    const float rate = (1.0f / supersamples);
	    const float r1 = sx * rate + rate / 2.0f;
	    const float r2 = sy * rate + rate / 2.0f;
	    // スクリーン上の位置
	    const Vec screen_position = 
	      screen_center + 
	      screen_x * ((r1 + x) / width - 0.5f) +
	      screen_y * ((r2 + y) / height- 0.5f);
	    // レイを飛ばす方向
	    const Vec dir = normalize(screen_position - camera_position);
	    Hitpoint null_hp;
	    accumulated_radiance = accumulated_radiance +  direct_radiance(Ray(camera_position, dir), &rnd, 0) / (supersamples * supersamples);
	    (*hdr.image_ptr(x, y)) = (*hdr.image_ptr(x, y)) + accumulated_radiance;					
	  }
	}
      }
    }

    std::cout << "done" << std::endl;
  }
	
  hdr.save("result0.hdr");

  return 0;
}
