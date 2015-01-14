#ifndef TASEP_HPP
#define TASEP_HPP

#include<vector>

using state_type = std::vector<double>;

class result_buffer {
public:
  std::vector<state_type>& plist;
  std::vector<double>& tlist;
  result_buffer(std::vector<state_type>& plist, std::vector<double>& tlist): plist(plist), tlist(tlist) {}
  void operator() (const state_type& p, double t);
  void write_to_file(const char* fname);
};

class tasep{
public:
  std::vector<double> rate;
  tasep(std::vector<double> rate) : rate(rate) {}
  void operator() (const state_type &p, state_type &dpdt, double t);
  double numeric_steady_state(result_buffer& observer, state_type& p, double t0=0, double t1=2000, double steady_threshold=1e-4);
  double numeric_steady_state(state_type& p, double t0=0, double t1=2000, double steady_threshold=1e-4);
};

#endif
