#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <string>
#include <sys/stat.h>
#include <chrono>

//------aliases------//
using rid_t = size_t;
using tid2refid_t = std::unordered_map<std::string, rid_t>;
using state_type = std::vector<double>;
using time_point = std::chrono::time_point<std::chrono::system_clock>;
using time_period = std::chrono::duration<double>;

//------const------//
const double EPSILON = 1e-7;
const unsigned SEED = 619048235;

//------string manipulation------//
bool startswith(const std::string &prefix, const std::string &input);
void string_split(const std::string &s, const std::string& delimiter, std::vector<std::string>& words);
void string_lstrip(const std::string& prefix, std::string &input);
void string_rstrip(const std::string& suffix, std::string &input);
inline bool fileExists(const std::string& filename);

template<typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& v)
{ 
  for (auto vi: v) out<<vi<<" "; 
  return out;
}

template<class vector_class>
std::ostream& print_vec(vector_class v, std::ostream& out=std::cout)
{
  for(auto vi:v) out<<vi<<" ";
  return out;
}

//------file utils------//
template<typename T>
void save_data_to_file(const std::string& fn, const std::vector<T>& data)
{
  std::ofstream ofile(fn);
  for (auto num: data)
    ofile<<num<<"\t";
  ofile<<std::endl;
  ofile.close();
}

template<typename T>
void append_data_to_file(const std::string& fn, const std::vector<T>& data)
{
  std::ofstream ofile;
  if (fileExists(fn))
    ofile.open(fn,std::ofstream::app);
  else
    ofile.open(fn, std::ofstream::out);
  for (auto num: data)
    ofile<<num<<"\t";
  ofile<<std::endl;
  ofile.close();
}

template<typename T>
void append_value_to_file(const std::string& fn, T val, const char sep = '\n')
{
  std::ofstream ofile;
  if (fileExists(fn))
    ofile.open(fn,std::ofstream::app);
  else
    ofile.open(fn, std::ofstream::out);
  ofile<<val<<sep;
  ofile.close();
}

// Function: fileExists
/**
    Check if a file exists
    @param[in] filename - the name of the file to check
    @return    true if the file exists, else false
*/
inline bool fileExists(const std::string& filename)
{
  struct stat buf;
  return (stat(filename.c_str(), &buf) != -1);
}
#endif
