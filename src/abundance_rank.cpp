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



#include <vector>
#include <string>
#include <numeric>
#include <fstream>
#include <cmath>

#include "reference_info_builder.hpp"
#include "ribomap_profiler.hpp"
#include "abundance_rank.hpp"

using namespace std;

abundance_rank::abundance_rank(const ribo_profile& rprofile, const transcript_info& tinfo)
{
  vector<rid_t> refID_vec(rprofile.get_expressed_transcript_ids());
  for (size_t t=0; t!=rprofile.number_of_transcripts(); ++t) {
    auto& p = rprofile.get_read_assignments(t);
    double rabd = rprofile.get_tot_count(t);
    // double rabd_sum = accumulate(p.begin(), p.end(), double(0));
    // if (std::fabs(rabd - rabd_sum)/rabd_original > 1e-3) {
    //   cout<<"accumulate abundance not equal tot_count! "<<rabd<<" "<<rprofile.get_tot_count(t)<<" "<<t<<endl;
    // }
    // if no ribosome footprint counts on any location > 1
    // or if no ribosome footprint counts on more than 2 loci > 0
    // don't consider ranking it
    if ((not have_spike(p)) or (not have_multiple_spikes(p))) continue;
    string tid = tinfo.get_tid(refID_vec[t]);
    double tabd = rprofile.get_tot_abundance(t)*rprofile.len(t);
    tlist.emplace_back(transcript_rank{tid, tabd, rabd, 0, 0, 0});
  } 
}

void abundance_rank::rank_by_tabd(int bin_num)
{
  sort(tlist.begin(), tlist.end(), 
       [](const transcript_rank & a, const transcript_rank & b) -> 
       bool { return a.tabd < b.tabd; }
       );
  size_t tcnt(tlist.size());
  long bin_size(tcnt/bin_num);
  for (size_t i=0; i!=tcnt; ++i) {
    auto& t=tlist[i];
    int bin_rank = i/bin_size;
    bin_rank = (bin_rank==bin_num) ? bin_num-1 : bin_rank;
    t.trank = bin_rank; 
  }
}

void abundance_rank::rank_by_rabd(int bin_num)
{
  sort(tlist.begin(), tlist.end(),
       [](const transcript_rank & a, const transcript_rank & b) ->
       bool { return a.rabd < b.rabd; }
       );
  size_t tcnt(tlist.size());
  long bin_size(tcnt/bin_num);
  for (size_t i=0; i!=tcnt; ++i) {
    auto& t=tlist[i];
    int bin_rank = i/bin_size;
    bin_rank = (bin_rank==bin_num) ? bin_num-1 : bin_rank;
    t.rrank = bin_rank;
  }
}

void abundance_rank::sort_by_rank_diff()
{
  for (size_t i=0; i!=tlist.size(); ++i) {
    auto& t=tlist[i];
    t.drank = t.trank-t.rrank;
  }
  sort(tlist.begin(), tlist.end(),
       [] (const transcript_rank& a, const transcript_rank& b) ->
       bool { return a.drank < b.drank; }
       );
}

void abundance_rank::get_rank(int bin_num)
{
  rank_by_tabd(bin_num);
  rank_by_rabd(bin_num);
  sort_by_rank_diff();
}

void abundance_rank::write_diff_list(const string& file_core,int percent_diff)
{
  string abd_fn(file_core+"_abundant.list"), scrc_fn(file_core+"_scarce.list");
  ofstream abd_file(abd_fn), scrc_file(scrc_fn);
  for (size_t i=0; i!=tlist.size(); ++i) {
    auto& t=tlist[i];
    if (t.drank < -percent_diff)
      abd_file<<t.tid<<" "<<t.tabd<<" "<<t.rabd<<" "<<t.trank<<" "<<t.rrank<<" "<<t.drank<<" "<<endl;
    else
      break;
  }
  for (size_t i=tlist.size()-1; i!=0; --i) {
    auto& t=tlist[i];
    if (t.drank > percent_diff)
      scrc_file<<t.tid<<" "<<t.tabd<<" "<<t.rabd<<" "<<t.trank<<" "<<t.rrank<<" "<<t.drank<<" "<<endl;
    else
      break;
  }
}
