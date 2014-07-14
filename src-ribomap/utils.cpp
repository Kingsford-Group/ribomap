#include "utils.hpp"
using namespace std;

void string_split(const string &s, const string& delimiter, vector<string>& words)
{
  size_t pre_pos = 0;
  size_t pos = 0;
  while ((pos = s.find(delimiter,pre_pos)) != string::npos) {
    words.push_back(s.substr(pre_pos, pos-pre_pos));
    pre_pos = pos+delimiter.size();
  }
  pos = s.size();
  words.push_back(s.substr(pre_pos, pos-pre_pos));
}

void string_lstrip(const string& prefix, string &input)
{
  if (input.find(prefix)==0)
    input.erase(0,prefix.size()); 
}

void string_rstrip(const string& suffix, string &input)
{
  int i=input.find(suffix);
  if (i!=0 && i!=string::npos)
    input = input.substr(0,i);
}
