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



#ifndef BAM_PARSER_HPP
#define BAM_PARSER_HPP

#include <map>
#include <unordered_map>
#include <set>
#include <vector>
#include <string>

using namespace std;

//------class forward declarations------//
class ribo_profile;
class transcript_info;

enum read_t { UNKNOWN, UTR5, FRAME0, FRAME1, FRAME2, UTR3 };

struct position{
  unsigned refID;
  int start;
  int stop;
  read_t type;
};

struct fp_record{
  int count;
  vector<position> al_loci;
  set<string> seqs;
  bool used;
};

// typedefs
using rd_rec_map_t = map<string, fp_record>;
using fp_list_t = vector<fp_record>;
using rlen2psite_t = map<int, int>;

// function declaration
bool get_expressed_alignments_from_bam(rd_rec_map_t& rd_rec, const char *fn, const ribo_profile& profiler, const string& cnt_sep, int lmin, int lmax, bool useSecondary, bool useRC);
// convert read range to P-site codon or codon range
read_t base_range_to_codon_range(const position& bp, int cds_begin, int cds_end, int phase, position& cp);
read_t base_range_to_middle_codon(const position& bp, int cds_begin, int cds_end, int offset, position& cp);
bool expressed_read_codon_ranges_from_bam(fp_list_t& fp_codon_list, const char *fn, const transcript_info& tinfo, const ribo_profile& profiler, const string& cnt_sep, int lmin, int lmax, bool useSecondary, bool useRC, int offset);
void alignment_regions_to_codon_ranges(const rd_rec_map_t& fp_rec, const transcript_info& tinfo, fp_list_t& fp_codon_list, int offset);
// get read type from read range -- no converting to codon
read_t read_type_from_range(const position& ibp, int cds_begin, int cds_end, position& obp);
read_t read_type_of_psite(const position& ibp, int cds_begin, int cds_end, int offset, position& obp);
// universal offset given
bool expressed_read_bases_from_bam(fp_list_t& fp_rec_out, const char* fn, const transcript_info& tinfo, const ribo_profile& profiler, const string& cnt_sep, int lmin, int lmax, bool useSecondary, bool useRC, int offset);
void assign_P_site(const rd_rec_map_t& fp_rec_in, const transcript_info& tinfo, fp_list_t& fp_rec_out, int offset);
// read length to offset file given
bool expressed_read_bases_from_bam(fp_list_t& fp_rec_out, const char* fn, const transcript_info& tinfo, const ribo_profile& profiler, const string& cnt_sep, int lmin, int lmax, bool useSecondary, bool useRC, const char* offset_fn);
rlen2psite_t get_readlen_psite_map(const char* offset_fn);
void assign_P_site(const rd_rec_map_t& fp_rec_in, const transcript_info& tinfo, fp_list_t& fp_rec_out, rlen2psite_t rl2p);
#endif


  
