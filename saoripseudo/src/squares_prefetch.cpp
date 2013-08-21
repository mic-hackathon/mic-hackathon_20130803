#include <ctime>
#include <sys/time.h> 
#include <iostream>
#include <cstdlib>
#include <math.h>
using namespace std;

#define MAX 3
#define MIN 0


double gettimeofday_sec()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec * 1e-6;
}

unsigned long xor128() {
  static unsigned long x=123456789, y=362436069, z=521288629, w=88675123;
  unsigned long t=(pow(x,(x<<11)));
  x=y; y=z; z=w;
  return ( w=pow((pow(w,(w>>19))),(pow(t,(t>>8)))) );
}

int main()
{
  const size_t num = 1000000;
  const float max = 3.0f;
  const float min = 0.0f;

  float* x = new float[num];
  float* y = new float[num];

  float a = 0.0f;
  float b = 0.0f;

  float sumxy = 0.0f;
  float sumx = 0.0f;
  float sumy = 0.0f;
  float sumx2 = 0.0f;

  double t1, t2;

  t1 = gettimeofday_sec();

#pragma omp parallel for 
  for(size_t i=0;i<num;i++){
    x[i] = static_cast<float>(i) * 0.1;
    y[i] = MIN + xor128() % (MAX - MIN) * 1.0 * static_cast<float>(i);
  }

#pragma prefetch x
#pragma prefetch y  
#pragma omp parallel reduction(+:sumxy, sumx, sumy, sumx2)
  {
  for(size_t i=0;i<num;i++){
    sumxy += x[i] * y[i];
    sumx  += x[i];
    sumy  += y[i];
    sumx2 += pow(x[i],2);
  }
  }

  a = (num * sumxy - sumx * sumy) / (num * sumx2 - pow(sumx, 2));
  b = (sumx2 * sumy - sumxy * sumx) / (num * sumx2 - pow(sumx, 2));

  t2 = gettimeofday_sec();

  cout << "a=" << a << endl;
  cout << "b=" << b << endl;
  cout << "time=" << t2 -t1 << endl;

  return 0;
  
}
