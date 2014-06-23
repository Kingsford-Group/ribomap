#ifndef GENCODE_PARSER_HPP
#define GENCODE_PARSER_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <cereal/access.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cereal/types/unordered_map.hpp>
#include <seqan/seq_io.h>

using namespace std;
using tid2refid_t = unordered_map<string,int>;

const int MAX_EXON_NUM = 100;
const unordered_map<string, string> codon2aa = {{"TTT", "F"}, {"TTC", "F"}, {"TTA", "L"}, {"TTG", "L"}, {"TCT", "S"}, {"TCC", "S"}, {"TCA", "S"}, {"TCG", "S"}, {"TAT", "Y"}, {"TAC", "Y"}, {"TAA", "U"}, {"TAG", "U"}, {"TGT", "C"}, {"TGC", "C"}, {"TGA", "U"}, {"TGG", "W"}, {"CTT", "L"}, {"CTC", "L"}, {"CTA", "L"}, {"CTG", "L"}, {"CCT", "P"}, {"CCC", "P"}, {"CCA", "P"}, {"CCG", "P"}, {"CAT", "H"}, {"CAC", "H"}, {"CAA", "Q"}, {"CAG", "Q"}, {"CGT", "R"}, {"CGC", "R"}, {"CGA", "R"}, {"CGG", "R"}, {"ATT","I"}, {"ATC", "I"}, {"ATA", "I"}, {"ATG", "M"}, {"ACT", "T"}, {"ACC", "T"}, {"ACA", "T"}, {"ACG", "T"}, {"AAT", "N"}, {"AAC", "N"}, {"AAA", "K"}, {"AAG", "K"}, {"AGT", "S"}, {"AGC", "S"}, {"AGA", "R"}, {"AGG", "R"}, {"GTT", "V"},{"GTC", "V"}, {"GTA", "V"}, {"GTG", "V"}, {"GCT", "A"}, {"GCC", "A"}, {"GCA", "A"}, {"GCG", "A"}, {"GAT", "D"}, {"GAC", "D"}, {"GAA", "E"}, {"GAG", "E"}, {"GGT", "G"}, {"GGC", "G"}, {"GGA", "G"},{"GGG", "G"}};

void encode_peptide(const string& tseq, string& pseq);

struct tprop{
  string tid;
  string gid;
  string tname;
  string gname;
  bool strand;
  int phase;
  int exon_num;
  int start;
  int stop;
  string chrm;
  int plen;
  tprop() {}
  tprop(const tprop& t): tid(t.tid), gid(t.gid), tname(t.tname), gname(t.gname), strand(t.strand), phase(t.phase), exon_num(t.exon_num), start(t.start), stop(t.stop), chrm(t.chrm), plen(t.plen) {}
  tprop(string tid, string gid, string tname, string gname, int start, int stop): tid(tid), gid(gid), tname(tname), gname(gname), strand(true), phase(0), exon_num(MAX_EXON_NUM), start(start), stop(stop), chrm(""), plen(0) {}
  template<class Archive>
  void serialize(Archive & archive) {
    archive(tid, gid, tname, gname, strand, phase, exon_num, start, stop, chrm, plen);
  }
  
  template<class Archive>
  static tprop* load_and_allocate( Archive & archive )
  {
    tprop tmp;
    archive(tmp.tid, tmp.gid, tmp.tname, tmp.gname, tmp.strand, tmp.phase, tmp.exon_num, tmp.start, tmp.stop, tmp.chrm, tmp.plen);
    return new tprop(tmp);
  }

};

class transcript_info{
public:
  vector<tprop> tlist;
  transcript_info(const char* tfname, const char* gtf_fname, const char* cereal_name = "../cache/tlist.cereal");
  size_t total_count() const { return tlist.size(); }
  uint32_t get_refID(const string& tid) const { return tid2refid.at(tid); }
  int cds_start(uint32_t refID) const { return tlist[refID].start-1; }
  int cds_stop(uint32_t refID) const { return tlist[refID].stop; }
  int cds_pep_len(uint32_t refID) const { return tlist[refID].plen;}
  string get_tid(uint32_t refID) const { return tlist[refID].tid; }
  int phase(uint32_t refID) const { return tlist[refID].phase; }
  template<class Archive>
  void serialize(Archive & archive) {
    archive(tlist, tid2refid);
  }
private:

  tid2refid_t tid2refid;
  bool get_info_from_fasta(const char* tfname);
  void build_tid_idx_map();
  bool get_info_from_gtf(const char* gtf_fname);
  void adjust_cds_ranges();
  void adjust_cds_ranges_pepfile(const char* pfname);
};

class fasta_reader{
public:
  seqan::FaiIndex faiIndex;
  fasta_reader(const char* fname);
  uint64_t length(uint32_t refID) const;
  string read_region(uint32_t refID, uint32_t start, uint32_t stop) const;
  string read_seq(uint32_t refID) const;
};
#endif
