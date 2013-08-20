#include <iostream>
#include <climits>
#include <cstring>
#include <cstdlib>
#include <sys/time.h>
#include <omp.h>

class Timer {
private:
  timeval t0_, t1_;
  std::string name_;
public:
  Timer(std::string name) : name_(name) {
    gettimeofday(&t0_, NULL);
  }

  double time_diff(timeval *begin, timeval *end) {
    return (double)(end->tv_sec - begin->tv_sec) + (end->tv_usec - begin->tv_usec) / 1000000.0;
  }
  
  double now() {
    gettimeofday(&t1_, NULL);
    return time_diff(&t0_, &t1_);
  }

  ~Timer() {
    //    std::cout << name_ << ": " << now() << " sec" << std::endl;
  }
};

class XorShift {
  unsigned int seed_[4];

  const float coefficient_;
 public:
  unsigned int next(void) { 
    const unsigned int t = seed_[0] ^ (seed_[0] << 11);
    seed_[0] = seed_[1]; 
    seed_[1] = seed_[2];
    seed_[2] = seed_[3];
    return seed_[3] = (seed_[3] ^ (seed_[3] >> 19)) ^ (t ^ (t >> 8)); 
  }

  float next01(void) {
    return (float)next() * coefficient_;
  }

  XorShift(const unsigned int initial_seed) :
    coefficient_(1.0f / (1.0f + UINT_MAX)){
    unsigned int s = initial_seed;
    for (int i = 1; i <= 4; i++){
      seed_[i-1] = s = 1812433253U * (s^(s>>30)) + i;
    }
  }
};
XorShift rnd(0);

enum Mode {
  ModeSeqWrite,
  ModeSeqRead,
  ModeSparseWrite,
  ModeSparseRead,
  ModeRandomWrite,
  ModeRandomRead,
};

void seq_write(unsigned char *data, const int length, const int num) {
  for (int i = 0; i < num; ++i) {
    const int idx = i % length;
    data[idx] = 0x12;
  }
}

unsigned char seq_read(unsigned char *data, const int length, const int num) {
  volatile unsigned char p = 0;
  for (int i = 0; i < num; ++i) {
    const int idx = i % length;
    p = data[idx];
  }
  return p;
}

void sparse_write(unsigned char *data, const int length, const int num) {
  for (int i = 0; i < num; ++i) {
    const int idx = (i + 65) % length;
    data[idx] = 0x12;
  }
}

unsigned char sparse_read(unsigned char *data, const int length, const int num) {
  volatile unsigned char p = 0;
  for (int i = 0; i < num; ++i) {
    const int idx = (i + 65) % length;
    p = data[idx];
  }
  return p;
}

void random_write(unsigned char *data, const int length, const int num) {
  for (int i = 0; i < num; ++i) {
    const int idx = rnd.next() % length;
    data[idx] = 0x12;
  }
}

unsigned char random_read(unsigned char *data, const int length, const int num) {
  volatile unsigned char p = 0;
  for (int i = 0; i < num; ++i) {
    const int idx = rnd.next() % length;
    p = data[idx];
  }
  return p;
}

int main(int argc, char **argv) {
  if (argc <= 3) {
    return 0;
  }
  const int num_thread = atoi(argv[1]);

  Mode mode = ModeRandomWrite;
  if (strcmp(argv[2], "seq_write") == 0)
    mode = ModeSeqWrite;
  else if (strcmp(argv[2], "seq_read") == 0)
    mode = ModeSeqRead;
  else if (strcmp(argv[2], "sparse_read") == 0)
    mode = ModeSparseRead;
  else if (strcmp(argv[2], "sparse_write") == 0)
    mode = ModeSparseWrite;
  else if (strcmp(argv[2], "random_write") == 0)
    mode = ModeRandomWrite;
  else if (strcmp(argv[2], "random_read") == 0)
    mode = ModeRandomRead;

  const int length = 1024 * 1024 * 1024; // 1024MB
  const int num = atoi(argv[3]);

  double proc_time[num_thread];
  unsigned char *data = new unsigned char[length];

  std::cout << "Thread Num: " << num_thread << std::endl;

  omp_set_num_threads(num_thread);
#pragma omp parallel for
  for (int i = 0; i < num_thread; ++i) {
    std::cout << "thread ID: " << i << std::endl;
    Timer timer("memory");
    switch (mode) {
    case ModeSeqWrite:
      seq_write(data, length, num);
      break;
    case ModeSeqRead:
      std::cout << seq_read(data, length, num) << std::endl;
      break;
      
    case ModeSparseWrite:
      sparse_write(data, length, num);
      break;
    case ModeSparseRead:
      std::cout << sparse_read(data, length, num) << std::endl;
      break;
      
    case ModeRandomWrite:
      random_write(data, length, num);
      break;
    case ModeRandomRead:
      std::cout << random_read(data, length, num) << std::endl;
      break;
    }
    proc_time[i] = timer.now();
  }

  for (int i = 0; i < num_thread; ++i) {
    std::cout << proc_time[i] << std::endl;
  }

  delete[] data;
  return 0;
}
