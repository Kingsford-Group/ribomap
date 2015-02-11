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



#ifndef ABUNDANCE_RANK_HPP
#define ABUNDANCE_RANK_HPP
#include <string>
#include <vector>

class ribo_profile;
class transcript_info;

struct transcript_rank {
  std::string tid;
  double tabd;
  double rabd;
  long trank;
  long rrank;
  long drank;
};

class abundance_rank {
public:
  const size_t size() const { return tlist.size(); }
  abundance_rank(const ribo_profile& rprofile, const transcript_info& tinfo);
  void get_rank(int bin_num);
  void write_diff_list(const string& file_core, int percent_diff);
private:
  std::vector<transcript_rank> tlist;
  void rank_by_tabd(int bin_num);
  void rank_by_rabd(int bin_num);
  void sort_by_rank_diff();
};

template<typename T>
bool have_spike(const std::vector<T>& p, double spike_level=1)
{
  for (auto pi: p) 
    if (pi>spike_level) return true;
  return false;
}

template<typename T>
bool have_multiple_spikes(const std::vector<T>& p, double spike_level=0)
{
  int i(0);
  for (auto pi: p)
    if (pi>spike_level) 
      ++i;
  if (i>2) return true;
  else return false;
}
#endif
