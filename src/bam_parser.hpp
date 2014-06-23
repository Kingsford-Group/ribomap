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

//------const------//
const int UTR3 = 1;
const int UTR5 = 2;
const int CDS = 0;

struct position{
  unsigned refID;
  int start;
  int stop;
  bool strand;
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
bool expressed_read_codon_ranges_from_bam(fp_list_t& fp_codon_list, const char *fn, const transcript_info& tinfo, const ribo_profile& profiler, int offset);
bool get_expressed_alignments_from_bam(rd_rec_map_t& rd_rec, const char *fn, const ribo_profile& profiler);
void alignment_regions_to_codon_ranges(const rd_rec_map_t& fp_rec, const transcript_info& tinfo, fp_list_t& fp_codon_list, int offset);
int base_range_to_middle_codon(const position& bp, int offset, int cds_begin, int cds_end, int phase, position& cp);
#endif


  
