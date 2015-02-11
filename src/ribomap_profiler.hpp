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



#ifndef RIBOMAP_PROFILER_HPP
#define RIBOMAP_PROFILER_HPP

#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include <iostream>

#include <boost/dynamic_bitset.hpp>
#include "bam_parser.hpp"
#include "utils.hpp"

using namespace std;

//------class forward declarations------//
class transcript_info;

//------classes------//
//------profile info per transcript------//
struct tprofile{
  // expected read count on each position 
  double tot_count;
  vector<double> count;
  // transcript abundance per base from sailfish
  double tot_abundance;
};

// // SHOULD GET RID OF
// //------expected count info per read per transcript-location pair------//
// struct read_loc_count {
//   rid_t tid;
//   rid_t loc;
//   double count;
// };
// //------count info per read------//
// struct read_count { 
//   double tot_count;
//   vector<read_loc_count> loc_count_list;
// };

/***
 * class ribo_profile:
 * store the profile from the model, the expected read assignment, and perform the EM
 */
class ribo_profile{
public:
  // vector<read_count> read_count_list; 
  unordered_map<rid_t, rid_t> refID2pID; //ref index from fasta --> profile index
  ribo_profile(const transcript_info& tinfo, const char* fname, const string& filetype, double abundance_cutoff = 1);
  bool assign_reads(fp_list_t& fp_base_list, const unordered_set<int>& type);
  vector<rid_t> get_expressed_transcript_ids() const;
  bool is_expressed(rid_t refID) const { return nonzero_abundance_vec[refID]; }
  const vector<double>& get_read_assignments(rid_t t) const { return profile[t].count;}
  double get_count(rid_t t, rid_t i) const { return profile[t].count[i]; }
  double get_tot_abundance(rid_t t) const { return profile[t].tot_abundance; }
  double get_tot_count(rid_t t) const { return profile[t].tot_count; }
  rid_t get_transcript_index(rid_t refID) const { return refID2pID.at(refID); }
  size_t len(rid_t t) const { return profile[t].count.size(); }
  size_t number_of_transcripts() const { return profile.size(); }
  bool check_tot_count(const fp_list_t& fp_base_list) const;
  // // Need to get rid of
  // bool initialize_read_count(const fp_list_t& fp_codon_list, bool normalize = true);
  // bool assign_reads();
  // void update_count_profile();
  // bool single_map_read_count(const fp_list_t& fp_codon_list);
private:
  vector<tprofile> profile;
  // bit vector: 1--transcript has non-zero abundance; 0--zero abundance
  boost::dynamic_bitset<> nonzero_abundance_vec; 
  double total_read_count;
  void sailfish_parser(const transcript_info& tinfo, const char* sf_fname, double abundance_cutoff);
  void cufflinks_parser(const transcript_info& tinfo, const char* cl_fname, double abundance_cutoff);
  void express_parser(const transcript_info& tinfo, const char* ep_fname, double abundance_cutoff);
  void build_read_prob(const transcript_info& tinfo);
  void reset_read_count();
  void include_abundant_transcript(rid_t refID) { nonzero_abundance_vec.set(refID); }
  void add_read_count(rid_t t, rid_t i, double count) { profile.at(t).count.at(i) += count; }
  void add_tot_count(rid_t t, double count) { profile[t].tot_count += count; }
};

#endif
