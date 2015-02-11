/**
	This file is part of the Ribomap suite
	-- an automatic pipeline for quantifying isoform-level
	ribosome profiles from ribosome profiling data.


	Copyright 2015 Hao Wang

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

	    http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
**/



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
