#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>
#include <utility>
#include <immintrin.h>


int solve0(std::vector<int>& p, std::vector<int>& w, int c)
{
  int n = p.size();
  std::vector<int> dp(c + 1, 0);

  for (int j = 0; j < n; ++j) {
    int pj = p[j];
    int wj = w[j];
    for (int k = c; k >= wj; --k) {
      int v0 = dp[k - wj] + pj;
      if (v0 > dp[k]) {
        dp[k] = v0;
      }
    }
  }
  return dp.back();
}


int solve1(std::vector<int>& p, std::vector<int>& w, int c)
{
  int n = p.size();
  std::vector<int> dp(c + 1);
  std::vector<int> dp0(c + 1, 0);

  for (int j = 0; j < n; ++j) {
    int pj = p[j];
    int wj = w[j];
    for (int k = 0; k < wj; ++k) {
      dp[k] = dp0[k];
    }
    for (int k = wj; k <= c; ++k) {
      dp[k] = std::max(dp0[k], dp0[k - wj] + pj);
    }
    dp.swap(dp0);
  }
  return dp0.back();
}


int solve2(std::vector<int>& p, std::vector<int>& w, int c)
{
  int n = p.size();
  std::vector<int> dp(c + 1);
  std::vector<int> dp0(c + 1, 0);

  for (int j = 0; j < n; ++j) {
    int pj = p[j];
    int wj = w[j];
#pragma omp parallel for
    for (int k = 0; k <= c; ++k) {
      int v;
      if (k < wj || dp0[k] >= (v = dp0[k - wj] + pj)) {
        dp[k] = dp0[k];
      } else {
        dp[k] = v;
      }
    }
    dp.swap(dp0);
  }
  return dp0.back();
}


int solve3(std::vector<int>& p, std::vector<int>& w, int c)
{
  int n = p.size();
  int vecsize = int(std::ceil(float(c + 1) / 16) * 16);
  int* dp = (int*)_mm_malloc(vecsize * sizeof(int), 512);
  int* dp0 = (int*)_mm_malloc(vecsize * sizeof(int), 512);
  int hoge[] __attribute__((aligned(64))) = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};

  for (int j = 0; j < n; ++j) {
    __m512i pj = _mm512_set1_epi32(p[j]);
    __m512i wj = _mm512_set1_epi32(w[j]);
    __m512i index = _mm512_set_epi32(15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    index = _mm512_add_epi32(index, wj);
    __m512i step = _mm512_set1_epi32(16);
    std::copy(dp0, dp0 + w[j], dp);
#pragma omp parallel for
    for (int k = w[j]; k <= c; k += 16) {
      //std::cout << "k = " << k << std::endl;
      //_mm512_store_epi32(hoge, index);
      //std::cout << "YamadaIsDead" << k << std::endl;
      //for (int i = 0; i < 16; ++i) {
      //  std::cout << hoge[i] << std::endl;
      //}
      __m512i v0 = _mm512_i32gather_epi32(index, dp0, 4);
      __m512i index2 = _mm512_sub_epi32(index, wj);
      __m512i v1 = _mm512_i32gather_epi32(index2, dp0, 4);
      v1 = _mm512_add_epi32(v1, pj);
      v1 = _mm512_max_epi32(v0, v1);
      //_mm512_store_epi32(dp+ k, v1);
      _mm512_i32scatter_epi32(dp, index, v1, 4);
      index = _mm512_add_epi32(index, step);
    }
    {
      int* tmp = dp;
      dp = dp0;
      dp0 = tmp;
    }
  }
  return dp0[c];
}


int solve(std::vector<int>& p, std::vector<int>& w, int c)
{
  return solve2(p, w, c);
}


int main(int argc, char* argv[])
{
  int n;
  int c;

  std::cin >> n >> c;

  std::vector<int> p(n);
  std::vector<int> w(n);

  for (int j = 0; j < n; ++j) {
    std::cin >> p[j] >> w[j];
  }

  int result = solve(p, w, c);

  int count = 1;
  auto start = std::chrono::system_clock::now();
  for (int i = 0; i < count; ++i) {
    solve(p, w, c);
  }
  auto stop = std::chrono::system_clock::now();
  auto time = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

  std::cout << "optimal value = " << result << std::endl;
  std::cout << "time = " << (double)time.count() / 1000 / count << std::endl;
  return 0;
}
