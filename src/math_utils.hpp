#ifndef MATH_UTILS_HPP
#define MATH_UTILS_HPP
#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>

template<typename T>
T median(const std::vector<T>& in_vec)
{
  std::vector<T> vec(in_vec);
  size_t size = vec.size();
  std::sort(vec.begin(), vec.end());
  return (size % 2 != 0) ? vec[size / 2] : 
    (vec[size / 2 - 1] + vec[size / 2]) / (T)2;
}

template<typename T>
void normalize_vec(std::vector<T>& num_vec)
{
  double vec_sum(0);
  for (auto num: num_vec)
    vec_sum += num;
  for (auto& num: num_vec)
    num /= vec_sum;
}

template<typename T>
void smooth_vec(std::vector<T>& num_vec, double epsilon=1e-5)
{
  for (auto& num: num_vec) num += epsilon;
}

template<typename T>
double euclidean_dist(const std::vector<T>& vec_a, const std::vector<T>& vec_b)
{
  double esum(0);
  assert(vec_a.size()==vec_b.size());
  for (size_t i=0; i!=vec_a.size(); ++i)
    esum += std::pow(vec_a[i]-vec_b[i], 2);
  return std::sqrt(esum);
}

template<typename T>
T l2norm(const std::vector<T>& nums)
{
  T sqr_sum = 0;
  for (auto x: nums)
    sqr_sum += std::pow(x,2);
  return std::sqrt(sqr_sum);
}

#endif
