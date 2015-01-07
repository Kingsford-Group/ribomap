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
};

// typedefs
using rd_rec_map_t = map<string, fp_record>;
using fp_list_t = vector<fp_record>;

// function declaration
bool get_expressed_alignments_from_bam(rd_rec_map_t& rd_rec, const char *fn, const ribo_profile& profiler, const string& cnt_sep);
// convert read range to P-site codon or codon range
bool expressed_read_codon_ranges_from_bam(fp_list_t& fp_codon_list, const char *fn, const transcript_info& tinfo, const ribo_profile& profiler, int offset, const string& cnt_sep);
void alignment_regions_to_codon_ranges(const rd_rec_map_t& fp_rec, const transcript_info& tinfo, fp_list_t& fp_codon_list, int offset);
read_t base_range_to_codon_range(const position& bp, int cds_begin, int cds_end, int phase, position& cp);
read_t base_range_to_middle_codon(const position& bp, int cds_begin, int cds_end, int offset, position& cp);
// get read type from read range -- no converting to codon
bool expressed_read_bases_from_bam(fp_list_t& fp_codon_list, const char* fn, const transcript_info& tinfo, const ribo_profile& profiler, int offset, const string& cnt_sep);
// universal offset given
void assign_P_site(const rd_rec_map_t& fp_rec_in, const transcript_info& tinfo, fp_list_t& fp_rec_out, int offset);
read_t read_type_from_range(const position& ibp, int cds_begin, int cds_end, position& obp);
read_t read_type_of_psite(const position& ibp, int cds_begin, int cds_end, int offset, position& obp);
#endif


  
