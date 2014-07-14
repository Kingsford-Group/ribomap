#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <vector>
#include <string>

void string_split(const std::string &s, const std::string& delimiter, std::vector<std::string>& words);
void string_lstrip(const std::string& prefix, std::string &input);
void string_rstrip(const std::string& suffix, std::string &input);

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v)
{ 
  for (auto num: v) out<<num<<" "; 
  return out;
}

#endif
